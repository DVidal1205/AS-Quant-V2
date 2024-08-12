import sys
import time
from methods.print_and_log import print_and_log
from methods.methods import SamtoTextParallel, SamtoTextSequential
import os
import glob
from collections import defaultdict
from multiprocessing import cpu_count
import csv


# Function to preprocess the GTF file
def preprocess_gene_ids(gtf, ref, outf, annotation_file, logger):
    # Step 0: Correct gene_id in GTF file using refFlat file
    print_and_log("----------------------------------------------------------------", logger)
    print_and_log("| STEP 0: Correcting GTF and creating annotation file...       |", logger)
    print_and_log("----------------------------------------------------------------", logger)

    # Track time
    start = time.monotonic()

    # Create a dictionary to store gene_id to gene_name mappings from refFlat file
    hash_dict = {}
    with open(ref, "r") as ref_file:
        for line in ref_file:
            info = line.split()
            # info[1] is gene_id and info[0] is gene_name in refFlat file
            hash_dict[info[1]] = info[0].replace('"', "").replace(";", "")

    # Dictionaries to store data for each transcript
    # transcripts (defaultdict)
    # |
    # |-- transcript_id_1 (str) -> {  "exons": [(start1, end1), (start2, end2), ...],
    # |                               "cds": [(start1, end1), (start2, end2), ...],
    # |                               "chrom": "chromosome",
    # |                               "strand": "+/-",
    # |                               "gene_id": "geneID_1",
    # |                               "gene_name": "geneName_1" }
    # |
    transcripts = defaultdict(lambda: {"exons": [], "cds": [], "chrom": None, "strand": None, "gene_name": None})

    # Open the GTF file for reading and the output file for writing
    count = 0
    drop_count = 0
    with open(gtf, "r") as gtf_file, open(outf, "w") as out_file:
        for line in gtf_file:
            line = line.strip()  # Remove leading/trailing whitespace
            info = line.split("\t")  # Split line by tab to get GTF fields

            # Check if the line contains gene_id and transcript_id
            if 'gene_id "' in info[8] and 'transcript_id "' in info[8]:
                gene_id = ""
                # Extract gene_id and transcript_id from the attributes field (info[8])
                g_id = info[8].split('gene_id "')[1].split('";')[0]
                g_id = g_id.split(".")[0]

                # Check for XM, XR, and YP prefixes in gene_id, and filter them out
                if g_id.startswith("XM") or g_id.startswith("XR") or g_id.startswith("YP"):
                    drop_count += 1
                    continue

                # If gene_id is in refFlat hash dictionary, use the corresponding gene_name
                if g_id in hash_dict:
                    gene_id = hash_dict[g_id]
                else:
                    # If gene_id is not in refFlat, drop the transcript
                    drop_count += 1
                    continue

                # Update the attributes field with the new or original gene_id
                info[8] = f'gene_id "{gene_id}"; transcript_id "{g_id}";'

                # Parse GTF line and store relevant information for annotation file
                chrom = info[0]
                feature_type = info[2]
                start = int(info[3])
                end = int(info[4])
                strand = info[6]
                transcript_id = g_id
                if feature_type == "exon":
                    transcripts[transcript_id]["exons"].append((start, end))
                elif feature_type == "CDS":
                    transcripts[transcript_id]["cds"].append((start, end))
                transcripts[transcript_id]["chrom"] = chrom
                transcripts[transcript_id]["strand"] = strand
                transcripts[transcript_id]["gene_name"] = gene_id

            else:
                print(f"Wrong line in GTF file:\n{line}", file=sys.stderr)
                sys.exit(1)

            # Write the updated line to the output file
            count += 1
            out_file.write("\t".join(info) + "\n")

    # Write to the annotation CSV
    with open(annotation_file, "w", newline="") as csv_file:
        writer = csv.writer(csv_file, delimiter="\t")

        # Write the header of the CSV file
        writer.writerow(["", "gene", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "txName"])

        # Write the data for each transcript to the CSV file
        for transcript_id, data in transcripts.items():
            # Grab all of the exon start and end positions
            exon_starts = ",".join(str(start) for start, end in data["exons"])
            exon_ends = ",".join(str(end) for start, end in data["exons"])
            # Get the minimum exon read from start and end positions
            tx_start = min(start for start, end in data["exons"])
            # Get the maximum exon read from start and end positions
            tx_end = max(end for start, end in data["exons"])
            # Get the minimum CDS read from start and end positions
            cds_start = min(start for start, end in data["cds"]) if data["cds"] else tx_start
            # Get the maximum CDS read from start and end positions
            cds_end = max(end for start, end in data["cds"]) if data["cds"] else tx_end
            writer.writerow([transcript_id, data["gene_name"], data["chrom"], data["strand"], tx_start, tx_end, cds_start, cds_end, len(data["exons"]), exon_starts, exon_ends, transcript_id])

    # Print the number of transcripts without a gene ID to stderr
    end = time.monotonic()
    print_and_log(f"{count} transcripts processed ", logger)
    print_and_log(f"{drop_count} transcripts dropped ", logger)



# Function to preprocess the BAM files using samtools
def preprocess_bam_files(input1_dir, input2_dir, species, samtools, cores, logger):
    # Step 1: Preprocess the BAM files using samtools
    print_and_log("----------------------------------------------------------------", logger)
    print_and_log("| STEP 2: Preprocessing BAM files using Samtools...            |", logger)
    print_and_log("----------------------------------------------------------------", logger)

    start = time.perf_counter()

    # First determine the chromosomes to use based off species
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
    if species == "hg38" or species == "hg19":
        chromosomes = chromosomes_h
    elif species == "mm10":
        chromosomes = chromosomes_m
    else:
        print_and_log("ERROR: SPECIES NOT FOUND. PLEASE SELECT: hg38, hg19, or mm10", logger)
        sys.exit()

    # Handle Cores
    parallel = True

    # Handle MAX and NULL flags
    if cores.upper() == "MAX" or cores.upper() == "NULL":
        if cores.upper() == "MAX":
            cores = cpu_count()
            print("Running AS-Quant with", cores, "cores...")
        else:
            cores = 1
            parallel = False

    # Grab the number of cores to use from the command line, if within the available range. Otherwise, exit.
    if int(cores) > cpu_count():
        cores = cpu_count()

    # If specified cores is equal to 1, then run the code in sequential mode
    if int(cores) == 1:
        parallel = False

    # Generate the coverage files for each chromosome using samtools
    current = os.getcwd()
    os.chdir(input1_dir)
    for file1 in glob.glob("*.bam"):
        if parallel:
            SamtoTextParallel(input1_dir, current, file1, chromosomes, samtools, cores, logger)
        else:
            SamtoTextSequential(input1_dir, current, file1, chromosomes, samtools, logger)
    os.chdir(os.path.join(current, input2_dir))
    for file2 in glob.glob("*.bam"):
        if parallel:
            SamtoTextParallel(input2_dir, current, file2, chromosomes, samtools, cores, logger)
        else:
            SamtoTextSequential(input2_dir, current, file2, chromosomes, samtools, logger)
    os.chdir(current)

    print_and_log(f"Step 2 Completed in {time.perf_counter() - start:.2f} seconds", logger)
    print_and_log("----------------------------------------------------------------", logger)
    logger.handlers[0].stream.write("\n\n")
    print()
