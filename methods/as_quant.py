import glob
import os
import warnings
import pandas as pd
from methods.methods import Count_pvalue_replicates, Count_pvalue, Find_Novel_splicing_events, Find_splicing_events, MakeFullDictionary, merge_ChromDict
from methods.print_and_log import print_and_log
import sys
from multiprocess import Pool, cpu_count
import time

warnings.simplefilter(action="ignore", category=FutureWarning)


def list_dirs(path):
    return [os.path.basename(x) for x in filter(os.path.isdir, glob.glob(os.path.join(path, "*")))]


def as_quant(species, output_dir, as_events, input1_dir, input2_dir, annotation_file, method, novel, cores, logger):
    # Step 3: Run AS-Quant
    print_and_log("----------------------------------------------------------------", logger)
    print_and_log("| STEP 3: Running AS-Quant...                                  |", logger)
    print_and_log("----------------------------------------------------------------", logger)

    start = time.perf_counter()

    chromosomes_h = [
        "chr1",
        "chr2",
        "chr3",
        "chr4",
        "chr5",
        "chr6",
        "chr7",
        "chr8",
        "chr9",
        "chr10",
        "chr11",
        "chr12",
        "chr13",
        "chr14",
        "chr15",
        "chr16",
        "chr17",
        "chr18",
        "chr19",
        "chr20",
        "chr21",
        "chr22",
        "chrX",
        "chrY",
    ]
    chromosomes_m = [
        "chr1",
        "chr2",
        "chr3",
        "chr4",
        "chr5",
        "chr6",
        "chr7",
        "chr8",
        "chr9",
        "chr10",
        "chr11",
        "chr12",
        "chr13",
        "chr14",
        "chr15",
        "chr16",
        "chr17",
        "chr18",
        "chr19",
        "chrX",
        "chrY",
    ]
    target_AS = ["SE", "RI", "MXE", "A3SS", "A5SS"]

    ########################################################################
    # Read in all of the command line arguments and set the default values #
    ########################################################################

    # Handle the method
    if method == "ranksum":
        count1 = len(glob.glob1(input1_dir, "*.bam"))
        count2 = len(glob.glob1(input2_dir, "*.bam"))
        if count1 < 2 or count2 < 2:
            print("Please provide multiple samples/replicates in each group to run ranksum test, otherwise select chisquare.")
            sys.exit()

    # Handle the cores
    parallel = True

    # Handle MAX and NULL flags
    if cores.upper() == "MAX" or cores.upper() == "NULL":
        if cores.upper() == "MAX":
            cores = cpu_count()
        else:
            cores = 1
            parallel = False

    # Grab the number of cores to use from the command line, if within the available range. Otherwise, exit.
    if int(cores) > cpu_count():
        cores = cpu_count()

    # If specified cores is equal to 1, then run the code in sequential mode
    if int(cores) == 1:
        parallel = False

    ########################################################################
    # Generate the coverage files for each chromosome based on species     #
    ########################################################################

    # Grab the group names from the input directories
    g1_name, g2_name = os.path.basename(input1_dir), os.path.basename(input2_dir)

    # Determine the chromosomes based on the species
    if species == "hg38" or species == "hg19":
        chromosomes = chromosomes_h
    elif species == "mm10":
        chromosomes = chromosomes_m
    else:
        print("Species not found. Please select among hg38, hg19 or mm10")
        sys.exit()

    # Load the annotation file.
    s1_namelist, s2_namelist = list_dirs(input1_dir), list_dirs(input2_dir)
    ann_df = pd.read_csv(annotation_file, delimiter="\t", index_col=0)

    # Convert the whole annotation into a dictionary for faster use
    print_and_log("Creating the chromosome dictionary...", logger)
    start = time.perf_counter()
    ChromDict = MakeFullDictionary(ann_df, chromosomes)
    print_and_log(f"Chromosome dictionary created in {time.perf_counter() - start:.2f} seconds", logger)

    # Merge the Exon Intervals
    print_and_log("Merging the Exon intervals...", logger)
    start = time.perf_counter()
    ChromDict_merged = merge_ChromDict(ChromDict, chromosomes)
    print_and_log(f"Exon intervals merged in {time.perf_counter() - start:.2f} seconds", logger)

    ########################################################################
    # Run the AS-Quant pipeline to detect significant splicing events      #
    ########################################################################

    print_and_log("Detecting Alternative Splicing Events...", logger)
    start = time.perf_counter()

    def execute_find_splicing_events(args):
        ChromDict_merged, chromosomes, AS, input_dir, species, sample, output_dir, as_events = args
        Find_splicing_events(ChromDict_merged, chromosomes, AS, input_dir, species, sample, output_dir, as_events)

    def parallel_find_splicing_events(ChromDict_merged, chromosomes, target_AS, input1_dir, input2_dir, species, s1_namelist, s2_namelist, output_dir, cores, as_events):
        tasks = []
        for AS in target_AS:
            as_start = time.perf_counter()
            print_and_log(f"Detecting {AS} events...", logger)
            for sample in s1_namelist:
                tasks.append((ChromDict_merged, chromosomes, AS, input1_dir, species, sample, output_dir, as_events))

            for sample in s2_namelist:
                tasks.append((ChromDict_merged, chromosomes, AS, input2_dir, species, sample, output_dir, as_events))

        with Pool(int(cores)) as pool:
            pool.map(execute_find_splicing_events, tasks)

    def execute_find_novel_splicing_events(args):
        (ChromDict_merged, ChromDict, chromosomes, AS, input_dir, species, sample, output_dir, as_events) = args
        Find_Novel_splicing_events(ChromDict_merged, ChromDict, chromosomes, AS, input_dir, species, sample, output_dir, as_events)

    def parallel_find_novel_splicing_events(
        ChromDict_merged,
        ChromDict,
        chromosomes,
        target_AS,
        input1_dir,
        input2_dir,
        species,
        s1_namelist,
        s2_namelist,
        output_dir,
        cores,
    ):
        tasks = []
        for AS in target_AS:
            as_start = time.perf_counter()
            print_and_log(f"Detecting {AS} events...", logger)

            for sample in s1_namelist:
                tasks.append((ChromDict_merged, ChromDict, chromosomes, AS, input1_dir, species, sample, output_dir, as_events))

            for sample in s2_namelist:
                tasks.append((ChromDict_merged, ChromDict, chromosomes, AS, input2_dir, species, sample, output_dir, as_events))

        with Pool(int(cores)) as pool:
            pool.map(execute_find_splicing_events, tasks)

        print_and_log(f"{AS} events found in {time.perf_counter() - as_start:.2f} seconds", logger)

    if novel.upper() == "YES":
        if parallel:
            parallel_find_novel_splicing_events(
                ChromDict_merged,
                ChromDict,
                chromosomes,
                target_AS,
                input1_dir,
                input2_dir,
                species,
                s1_namelist,
                s2_namelist,
                output_dir,
                cores,
            )
        else:
            # OLD: Sequential version of Find_Novel_splicing_events
            target_AS = ["All"]
            for AS in target_AS:
                for sample in s1_namelist:
                    Find_Novel_splicing_events(
                        ChromDict_merged,
                        ChromDict,
                        chromosomes,
                        AS,
                        input1_dir,
                        species,
                        sample,
                        output_dir,
                    )
                for sample in s2_namelist:
                    Find_Novel_splicing_events(
                        ChromDict_merged,
                        ChromDict,
                        chromosomes,
                        AS,
                        input2_dir,
                        species,
                        sample,
                        output_dir,
                    )

    else:
        if parallel:
            parallel_find_splicing_events(ChromDict_merged, chromosomes, target_AS, input1_dir, input2_dir, species, s1_namelist, s2_namelist, output_dir, cores, as_events)
        else:
            # OLD: Sequential version of Find_splicing_events
            for AS in target_AS:
                as_start = time.perf_counter()
                print_and_log(f"Detecting {AS} events...", logger)
                for sample in s1_namelist:
                    Find_splicing_events(ChromDict_merged, chromosomes, AS, input1_dir, species, sample, output_dir, as_events)
                for sample in s2_namelist:
                    Find_splicing_events(ChromDict_merged, chromosomes, AS, input2_dir, species, sample, output_dir, as_events)
                print_and_log(f"{AS} events found in {time.perf_counter() - start:.2f} seconds", logger)

    print_and_log(f"Found alternative splicing events in {time.perf_counter() - start:.2f} seconds", logger)
    print_and_log("Writing to output...", logger)

    ########################################################################
    # Write the output to an Excel file 								   #
    ########################################################################

    writer_out = pd.ExcelWriter(
        os.path.join(output_dir, "asquant_" + g1_name + "_Vs_" + g2_name + ".xlsx"),
        engine="xlsxwriter",
    )
    for AS in target_AS:
        if method.lower() == "ranksum":
            Count_pvalue_replicates(AS, output_dir, s1_namelist, s2_namelist)
        else:
            Count_pvalue(AS, output_dir, s1_namelist, s2_namelist, g1_name, g2_name)

        df = pd.read_csv(os.path.join(output_dir, AS + "_Output.csv"), delimiter="\t")
        df.sort_values(by=["P_value"], ascending=True, inplace=True)
        df.to_excel(writer_out, sheet_name=AS, index=None, header=True)
        os.remove(os.path.join(output_dir, AS + "_Output.csv"))
    writer_out.save()

    print_and_log(f"Step 3 Completed in {time.perf_counter() - start:.2f} seconds", logger)
    print_and_log("----------------------------------------------------------------", logger)
    logger.handlers[0].stream.write("\n\n")
    print()
