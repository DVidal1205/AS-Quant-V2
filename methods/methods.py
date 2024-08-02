from operator import itemgetter
import pandas as pd
import time
from bisect import bisect_left
import os
from scipy.stats import chisquare
import numpy as np
from scipy import stats
import sys
from multiprocessing import Pool
from methods.print_and_log import print_and_log


class Stack:
    def __init__(self):
        self.items = []

    def size(self):
        return len(self.items)

    def isEmpty(self):
        return self.items == []

    def push(self, val):
        self.items.append(val)

    def top(self):
        if self.isEmpty():
            return None
        else:
            return self.items[self.size() - 1]

    def pop(self):
        if self.isEmpty():
            return None
        else:
            return self.items.pop()


def bi_contains(lst, item):
    return bisect_left(lst, item)


def MergeIntervals(inputlist):
    n = len(inputlist)
    inputlist.sort(key=itemgetter(1), reverse=False)

    st = Stack()
    st.push(inputlist[0])

    for i in range(1, n):
        stacktop = st.top()
        if inputlist[i][0] <= stacktop[1]:
            st.pop()
            st_st = stacktop[0]
            st_en = max(stacktop[1], inputlist[i][1])
            st.push((st_st, st_en))
        else:
            st.push(inputlist[i])

    mergedExList = []
    while True:
        if st.size() == 0:
            break
        stacktop = st.top()
        mergedExList.append(stacktop)
        st.pop()

    return mergedExList


def CountTotalReadCount(chrom, exList, bam_list, position_row):
    totalCount = 0
    for p in range(len(exList)):
        start = int(exList[p][0])
        end = int(exList[p][1])

        pos1 = bi_contains(position_row, start)
        pos2 = bi_contains(position_row, end)

        if pos1 < len(position_row) and pos2 < len(position_row):
            if int(bam_list[pos2][0]) != end:
                pos2 = pos2 - 1

            for t in range(pos1, pos2 + 1):
                read = int(bam_list[t][1])
                totalCount += read

    return totalCount


def writeResult(chrom, gene, start, end, bam_list, position_row, RC, mergedExListLength, writer_list):
    targetRC = CountTotalReadCount(chrom, [(start, end)], bam_list, position_row)
    targetLength = end - start + 1

    #### Avoiding divide by zero error #####
    if targetLength == 0:
        averageTargetRC = 0
    else:
        averageTargetRC = targetRC / targetLength

    if mergedExListLength == targetLength:
        averageRCothers = 0
    else:
        averageRCothers = (RC - targetRC) / (mergedExListLength - targetLength)

    writer_list.append(
        (
            chrom,
            gene,
            start,
            end,
            targetRC,
            targetLength,
            RC,
            mergedExListLength,
            RC - targetRC,
            mergedExListLength - targetLength,
            averageTargetRC,
            averageRCothers,
        )
    )

    return writer_list


def Find_Novel_splicing_events(ChromDict_merged, ChromDict, chromosomes, AS, input_dir, species_folder, sample, output_dir, as_events_dir):
    tt = time.time()
    AS_flag = []
    as_df = pd.read_csv(os.path.join(as_events_dir, species_folder, AS + ".csv"), delimiter="\t")
    writer_list = []
    output_columns = [
        "chrom",
        "geneName",
        "splicedExonStart",
        "splicedExonEnd",
        "splicedExonReadCount(rc)",
        "splicedExonLength(sl)",
        "othersExonsReadCount(RC)",
        "othersExonsLength(L)",
        "RC - rc",
        "L - sl",
        "splicedExonAverageReadCoverage(n)",
        "otherExonsAverageReadCoverage(N)",
    ]

    for chrom in chromosomes:
        # print("Starting:",chrom)
        tts = time.time()
        GeneDict = ChromDict[chrom]
        GeneDict_merged = ChromDict_merged[chrom]
        if os.path.getsize(os.path.join(input_dir, sample, chrom + ".txt")) > 0:
            bam_df = pd.read_csv(os.path.join(input_dir, sample, chrom + ".txt"), delimiter="\t")
            position_row = bam_df.iloc[:, 0].tolist()
            bam_list = bam_df.values.tolist()

            for gene in GeneDict.keys():
                exonList = list(set(GeneDict[gene.upper()]))
                if gene.upper() in GeneDict_merged:
                    mergedExList = GeneDict_merged[gene.upper()]

                mergedExListLength = 0
                for p in range(len(mergedExList)):
                    mergedExListLength += mergedExList[p][1] - mergedExList[p][0] + 1
                RC = CountTotalReadCount(chrom, mergedExList, bam_list, position_row)

                for ex in range(len(exonList)):
                    start, end = int(exonList[ex][0]), int(exonList[ex][1])

                    if (chrom, gene, start, end) not in AS_flag:
                        writer_list = writeResult(
                            chrom,
                            gene,
                            start,
                            end,
                            bam_list,
                            position_row,
                            RC,
                            mergedExListLength,
                            writer_list,
                        )
                        AS_flag.append((chrom, gene, start, end))

    df_out = pd.DataFrame(writer_list, columns=output_columns)
    df_out.to_csv(os.path.join(output_dir, sample + "_" + AS + ".csv"), sep="\t", index=False)


def Find_splicing_events(ChromDict_merged, chromosomes, AS, input_dir, species, sample, output_dir, as_events_dir):
    tt = time.time()
    AS_flag = []
    as_df = pd.read_csv(os.path.join(as_events_dir, species, AS + ".csv"), delimiter="\t")

    writer_list = []
    output_columns = [
        "chrom",
        "geneName",
        "splicedExonStart",
        "splicedExonEnd",
        "splicedExonReadCount(rc)",
        "splicedExonLength(sl)",
        "othersExonsReadCount(RC)",
        "othersExonsLength(L)",
        "RC - rc",
        "L - sl",
        "splicedExonAverageReadCoverage(n)",
        "otherExonsAverageReadCoverage(N)",
    ]

    for chrom in chromosomes:
        # print("Starting:",chrom)
        tts = time.time()
        GeneDict = ChromDict_merged[chrom]
        if os.path.getsize(os.path.join(input_dir, sample, chrom + ".txt")) > 0:
            bam_df = pd.read_csv(os.path.join(input_dir, sample, chrom + ".txt"), delimiter="\t")
            position_row = bam_df.iloc[:, 0].tolist()
            bam_list = bam_df.values.tolist()

            as_chr_rows = as_df[as_df["chrom"] == chrom]
            for ind1, t_row in as_chr_rows.iterrows():
                gene = t_row["gene"].strip().upper()
                if gene in GeneDict:
                    mergedExList = GeneDict[gene]

                mergedExListLength = 0
                for p in range(len(mergedExList)):
                    mergedExListLength += mergedExList[p][1] - mergedExList[p][0] + 1

                RC = CountTotalReadCount(chrom, mergedExList, bam_list, position_row)

                if AS in ["SE", "RI"]:
                    exonStart, exonEnd = t_row["exonStart"], t_row["exonEnd"]
                    if (chrom, gene, exonStart, exonEnd) not in AS_flag:
                        writer_list = writeResult(
                            chrom,
                            gene,
                            exonStart,
                            exonEnd,
                            bam_list,
                            position_row,
                            RC,
                            mergedExListLength,
                            writer_list,
                        )
                        AS_flag.append((chrom, gene, exonStart, exonEnd))

                elif AS == "MXE":
                    exon1Start, exon1End = t_row["exon1Start"], t_row["exon1End"]
                    exon2Start, exon2End = t_row["exon2Start"], t_row["exon2End"]
                    if (chrom, gene, exon1Start, exon1End) not in AS_flag:
                        writer_list = writeResult(
                            chrom,
                            gene,
                            exon1Start,
                            exon1End,
                            bam_list,
                            position_row,
                            RC,
                            mergedExListLength,
                            writer_list,
                        )
                        AS_flag.append((chrom, gene, exon1Start, exon1End))

                    if (chrom, gene, exon2Start, exon2End) not in AS_flag:
                        writer_list = writeResult(
                            chrom,
                            gene,
                            exon2Start,
                            exon2End,
                            bam_list,
                            position_row,
                            RC,
                            mergedExListLength,
                            writer_list,
                        )
                        AS_flag.append((chrom, gene, exon2Start, exon2End))

                else:
                    longExonStart, longExonEnd, shortExonStart, shortExonEnd, strand = (
                        t_row["longExonStart"],
                        t_row["longExonEnd"],
                        t_row["shortExonStart"],
                        t_row["shortExonEnd"],
                        t_row["strand"],
                    )
                    start, end = 0, 0
                    if AS == "A5SS":
                        if strand == "+":
                            start, end = longExonEnd + 1, shortExonEnd
                        else:
                            start, end = shortExonStart, longExonStart - 1

                    elif AS == "A3SS":
                        if strand == "+":
                            start, end = shortExonStart, longExonStart - 1
                        else:
                            start, end = longExonEnd + 1, shortExonEnd

                    if (chrom, gene, start, end) not in AS_flag:
                        writer_list = writeResult(
                            chrom,
                            gene,
                            start,
                            end,
                            bam_list,
                            position_row,
                            RC,
                            mergedExListLength,
                            writer_list,
                        )
                        AS_flag.append((chrom, gene, start, end))

    df_out = pd.DataFrame(writer_list, columns=output_columns)
    df_out.to_csv(os.path.join(output_dir, sample + "_" + AS + ".csv"), sep="\t", index=False)


def MakeFullDictionary(ann_df, chromosomes):
    ChromDict = {}
    for chrom in chromosomes:
        GeneDict = {}
        chr_rows = ann_df[ann_df["chrom"] == chrom]
        gene_list = list(set(chr_rows["gene"]))
        for gene in gene_list:
            gene_rows = chr_rows[chr_rows["gene"] == gene]
            exList = []
            for index, row in gene_rows.iterrows():
                exonCount = row["exonCount"]
                exonStarts = list(filter(None, row["exonStarts"].split(",")))
                exonEnds = list(filter(None, row["exonEnds"].split(",")))
                for i in range(exonCount):
                    st, en = int(exonStarts[i]), int(exonEnds[i])
                    if (st, en) not in exList:
                        exList.append((st, en))

            GeneDict[gene.strip().upper()] = exList

        ChromDict[chrom] = GeneDict

    return ChromDict


def merge_ChromDict(ChromDict, chromosomes):
    ChromDict_merged = {}
    for chrom in chromosomes:
        GeneDict_merged = {}
        GeneDict = ChromDict[chrom]
        for gene in GeneDict.keys():
            exonList = GeneDict[gene.upper()]
            mergedExonList = MergeIntervals(exonList)
            GeneDict_merged[gene.upper()] = mergedExonList
        ChromDict_merged[chrom] = GeneDict_merged

    return ChromDict_merged


def Count_pvalue_replicates(AS, output_dir, s1_namelist, s2_namelist, g1_name, g2_name):
    len_S1 = len(s1_namelist)
    len_S2 = len(s2_namelist)

    S1_list, S2_list = {}, {}
    writer_list = []
    output_columns = [
        "chrom",
        "geneName",
        "splicedExonStart",
        "splicedExonEnd",
        "P_value",
        "ChromRegionLong",
    ]

    for k, sample1 in enumerate(s1_namelist):
        S1_df = pd.read_csv(os.path.join(output_dir, sample1 + "_" + AS + ".csv"), delimiter="\t")
        S1_list[k] = list(S1_df.values.tolist())

    for p, sample2 in enumerate(s2_namelist):
        S2_df = pd.read_csv(os.path.join(output_dir, sample2 + "_" + AS + ".csv"), delimiter="\t")
        S2_list[p] = list(S2_df.values.tolist())
        AS_file_length = len(S2_df)

    #### need to check this part, whether it can read the values
    for i in range(1, AS_file_length):
        X, Y = [], []
        N1s, N2s = [], []
        for k in range(len_S1):
            full_list = S1_list[k]
            chrom = full_list[i][0]
            gene = full_list[i][1]
            start = full_list[i][2]
            end = full_list[i][3]

            n1 = float(full_list[i][10])
            N1 = float(full_list[i][11])
            N1s.append(N1)
            N1 = N1 + n1
            if N1 > 0:
                r1 = n1 / N1
                X.append(r1)

        for p in range(len_S2):
            full_list = S2_list[p]
            n2 = float(full_list[i][10])
            N2 = float(full_list[i][11])
            N2s.append(N2)
            N1 = N2 + n2
            if N2 > 0:
                r2 = n2 / N2
                Y.append(r2)

        ##### added N filter: 11.14.2022
        ##### both N should be >= 10
        if (min(np.mean(N1s), np.mean(N2s)) >= 10) and sum(np.array(X)) != 0 and sum(np.array(Y)) != 0:
            # stat, p_val = ranksums(X, Y)
            stat, p_val = stats.ttest_ind(X, Y)
            chrom_region_long = chrom + ":" + gene + ":" + str(start) + "-" + str(end)
            writer_list.append((chrom, gene, start, end, p_val, chrom_region_long))

    df_out = pd.DataFrame(writer_list, columns=output_columns)
    df_out.to_csv(os.path.join(output_dir, AS + "_Output.csv"), sep="\t", index=False)

    for sample1 in s1_namelist:
        os.remove(os.path.join(output_dir, sample1 + "_" + AS + ".csv"))
    for sample2 in s2_namelist:
        os.remove(os.path.join(output_dir, sample2 + "_" + AS + ".csv"))

    return


def Count_pvalue(AS, output_dir, s1_namelist, s2_namelist, g1_name, g2_name):
    len_S1, len_S2 = len(s1_namelist), len(s2_namelist)
    S1_list, S2_list = {}, {}

    writer_list = []
    output_columns = [
        "chrom",
        "geneName",
        "splicedExonStart",
        "splicedExonEnd",
        "P_value",
        "ratioDifference",
        "absoluteRatioDifference",
        g1_name + ":n1",
        g1_name + ":N1",
        g2_name + ":n2",
        g2_name + ":N2",
        "chromRegionLLong",
    ]

    for k, sample1 in enumerate(s1_namelist):
        S1_df = pd.read_csv(os.path.join(output_dir, sample1 + "_" + AS + ".csv"), delimiter="\t")
        S1_list[k] = list(S1_df.values.tolist())

    for p, sample2 in enumerate(s2_namelist):
        S2_df = pd.read_csv(os.path.join(output_dir, sample2 + "_" + AS + ".csv"), delimiter="\t")
        S2_list[p] = list(S2_df.values.tolist())
        AS_file_length = len(S2_df)

    for i in range(1, AS_file_length):
        S1_n_total, S1_N_total, S2_n_total, S2_N_total = 0.0, 0.0, 0.0, 0.0
        for k in range(len_S1):
            full_list = S1_list[k]
            S1_n_total += float(full_list[i][10])
            S1_N_total += float(full_list[i][11])

            chrom = full_list[i][0]
            gene = full_list[i][1]
            start = full_list[i][2]
            end = full_list[i][3]

        n1 = S1_n_total / len_S1
        N1 = S1_N_total / len_S1

        for p in range(len_S2):
            full_list = S2_list[p]
            S2_n_total += float(full_list[i][10])
            S2_N_total += float(full_list[i][11])

        n2 = S2_n_total / len_S2
        N2 = S2_N_total / len_S2

        ##### added N filter: 11.14.2022
        ##### both N should be >= 10 and no negative values of n and N are allowed
        if (min(N1, N2) >= 10) and (n1 < N1) and (n2 < N2):
            N1 = N1 + n1
            N2 = N2 + n2
            if N1 > 0 and N2 > 0:
                ratio_diff = (n1 / N1) - (n2 / N2)
                abs_ratio_diff = abs(ratio_diff)

                P0 = (n1 + n2) / (N1 + N2)
                n10 = N1 * P0
                n20 = N2 * P0
                exp = [n10, N1 - n10, n20, N2 - n20]
                if 0 not in exp:
                    res = chisquare([n1, N1 - n1, n2, N2 - n2], f_exp=exp, ddof=1)
                    chrom_region_long = chrom + ":" + gene + ":" + str(start) + "-" + str(end)
                    writer_list.append(
                        (
                            chrom,
                            gene,
                            start,
                            end,
                            res[1],
                            ratio_diff,
                            abs_ratio_diff,
                            n1,
                            N1 - n1,
                            n2,
                            N2 - n2,
                            chrom_region_long,
                        )
                    )

    df_out = pd.DataFrame(writer_list, columns=output_columns)
    df_out.to_csv(os.path.join(output_dir, AS + "_Output.csv"), sep="\t", index=False)

    for sample1 in s1_namelist:
        os.remove(os.path.join(output_dir, sample1 + "_" + AS + ".csv"))
    for sample2 in s2_namelist:
        os.remove(os.path.join(output_dir, sample2 + "_" + AS + ".csv"))

    return


# Function to process each chromosome that will be mapped to each core
def process_chromosome(args):
    samtools_dir, input_dir, current, bamfile_name, chrom, output_dir = args
    try:
        cmd2 = f"{samtools_dir} view -b {os.path.join(current, input_dir, bamfile_name)} {chrom} -o {os.path.join(current, output_dir, chrom + '.bam')}"
        cmd3 = f"{samtools_dir} pileup {os.path.join(current, output_dir, chrom + '.bam')} | cut -f 2,4 > {os.path.join(current, output_dir, chrom + '.txt')}"
        command = f"{cmd2}; {cmd3}"
        os.system(command)
    except ValueError:
        print(f"Read coverage file for chromosome {chrom} could not be generated")
        sys.exit()


# Function to generate read coverage files for each chromosome in parallel
def SamtoTextParallel(input_dir, current, bamfile_name, chromosomes, samtools_dir, cores, logger):
    start = time.perf_counter()

    output_dir = os.path.join(input_dir, os.path.splitext(bamfile_name)[0])
    last_dir = os.path.basename(os.path.normpath(bamfile_name)).split(".bam")[0]
    os.makedirs(last_dir, exist_ok=True)

    try:
        cmd1 = f"{samtools_dir} index {os.path.join(current, input_dir, bamfile_name)}"
        os.system(cmd1)
    except ValueError:
        print_and_log("Index file could not be generated", logger)
        sys.exit()

    args_list = [(samtools_dir, input_dir, current, bamfile_name, chrom, output_dir) for chrom in chromosomes]

    with Pool(int(cores)) as pool:
        pool.map(process_chromosome, args_list)

    for chrom in chromosomes:
        bam_path = os.path.join(current, output_dir, chrom + ".bam")
        if os.path.exists(bam_path):
            os.remove(bam_path)

    message = "Read coverage files generated for " + bamfile_name + " in " + f"{time.perf_counter() - start:.2f} seconds"
    print_and_log(message, logger)
    return


# Function to generate read coverage files for each chromosome sequentially
def SamtoTextSequential(input_dir, current, bamfile_name, chromosomes, samtools_dir, logger):
    start = time.perf_counter()

    output_dir = os.path.join(input_dir, os.path.splitext(bamfile_name)[0])
    last_dir = os.path.basename(os.path.normpath(bamfile_name)).split(".bam")[0]
    os.makedirs(last_dir, exist_ok=True)

    try:
        cmd1 = samtools_dir + " index " + os.path.join(current, input_dir, bamfile_name)
        os.system(cmd1)
    except ValueError:
        print_and_log("Index file could not be generated", logger)
        sys.exit()

    for chrom in chromosomes:
        cmd2 = samtools_dir + " view -b " + os.path.join(current, input_dir, bamfile_name) + " " + chrom + " -o " + os.path.join(current, output_dir, chrom + ".bam")
        cmd3 = samtools_dir + " pileup " + os.path.join(current, output_dir, chrom + ".bam") + " | cut -f 2,4 > " + os.path.join(current, output_dir, chrom + ".txt")  ### Need to use pileup, not mpileup
        command = cmd2 + ";" + cmd3
        try:
            os.system(command)
            os.system("rm " + os.path.join(current, output_dir, chrom + ".bam"))
        except ValueError:
            print_and_log("Read coverage file could not be generated", logger)
            sys.exit()

    message = "Read coverage files generated for " + bamfile_name + " in " + f"{time.perf_counter() - start:.2f} seconds"
    print_and_log(message, logger)
    return
