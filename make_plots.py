import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
from bisect import bisect_left
import sys
import os
import glob
from matplotlib import rcParams
from methods.methods import MakeFullDictionary, merge_ChromDict
import configparser


plt.rc("legend", **{"fontsize": 8})

rcParams.update(
    {
        "font.family": "arial",
    }
)

y_limit = 0


def bi_contains(lst, item):
    return bisect_left(lst, item)


def list_dirs(path):
    return [os.path.basename(x) for x in filter(os.path.isdir, glob.glob(os.path.join(path, "*")))]


def Generate_read_coverate_plot(ax, pathin, sample, chrom, geneID, start, end, startAll, endAll, group_name, number):
    bam_df = pd.read_csv(os.path.join(pathin, sample, chrom + ".txt"), delimiter="\t")
    position_row = bam_df.iloc[:, 0].tolist()
    bam_list = bam_df.values.tolist()

    ax = ax or plt.gca()
    x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    x_formatter.set_scientific(False)
    ax.xaxis.set_major_formatter(x_formatter)

    pos1 = bi_contains(position_row, startAll)
    pos2 = bi_contains(position_row, endAll)
    if int(bam_list[pos2][0]) != endAll:
        pos2 = pos2 - 1

    p = []
    c = []
    read = 0
    length = endAll - startAll + 1

    for t in range(length):
        p.append(t + startAll)
        c.append(0)

    for t in range(pos1, pos2 + 1):
        position = int(bam_list[t][0])
        read = int(bam_list[t][1])
        index = p.index(position)
        c[index] = read

    p = np.array(p)
    c = np.array(c)

    pos3 = bi_contains(p, start)
    pos4 = bi_contains(p, end)

    global y_limit
    m = max(c)
    if m > y_limit:
        y_limit = m

    if number == 1:
        caption = ax.fill_between(p, c, color="midnightblue", alpha=0.9, label=group_name)
    else:
        caption = ax.fill_between(p, c, color="crimson", alpha=0.9, label=group_name)

    ax.legend(handles=[caption])
    ax.fill_between(p[pos3 : pos4 + 1], c[pos3 : pos4 + 1], color="orange", alpha=0.9)
    ax.spines["top"].set_color("none")
    ax.spines["right"].set_color("none")
    ax.set_xlim(startAll, endAll)
    ax.autoscale(enable=True)
    ax.set_xticklabels([])
    ax.tick_params(axis="both", bottom=False, which="major", labelsize=8)

    return y_limit


def Generate_annotation_plot(ax, isoforms, exonCountList, exonStartList, exonEndList, start, end, startAll, endAll):
    ax = ax or plt.gca()
    x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    x_formatter.set_scientific(False)
    ax.xaxis.set_major_formatter(x_formatter)
    print("Total number of isoforms are: ", isoforms)

    ystart = 0
    height = 3
    for i in range(isoforms):
        if i >= 15:
            print("15 isoforms of this Gene is plotted.")
            break
        else:
            ax.hlines(
                y=(ystart + ystart + height) / 2,
                xmin=startAll,
                xmax=endAll,
                linewidth=1,
                color="skyblue",
                linestyle="--",
            )

            ecount = int(exonCountList[i])
            stList = exonStartList[i]
            enList = exonEndList[i]
            for p in range(ecount):
                ex_s = int(stList[p])
                ex_e = int(enList[p])
                width = int(enList[p]) - int(stList[p]) + 1

                rect = patches.Rectangle((ex_s, ystart), width, height, color="skyblue", fill=True)
                ax.add_patch(rect)

                if start >= ex_s and end <= ex_e:
                    width = end - start + 1
                    rect1 = patches.Rectangle((start, ystart), width, height, color="orange", alpha=0.9, fill=True)
                    ax.add_patch(rect1)
                elif start <= ex_s and end >= ex_e:
                    width = ex_e - ex_s + 1
                    rect1 = patches.Rectangle((ex_s, ystart), width, height, color="orange", alpha=0.9, fill=True)
                    ax.add_patch(rect1)

            ystart += 5

    ax.set_xlim(startAll, endAll)
    ax.autoscale(enable=True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.set_yticklabels([])
    ax.tick_params(left=False, axis="both", which="major", labelsize=10)

    return


def Process_user_inputs(
    region,
    ChromDict_merged,
    input1_dir,
    input2_dir,
    s1_namelist,
    s2_namelist,
    g1_name,
    g2_name,
    output_dir,
    ann_df,
):
    chrom, geneID, rng = region.split(":")
    start, end = rng.split("-")

    GeneDict_merged = ChromDict_merged[chrom]
    exList = GeneDict_merged[geneID.strip().upper()]

    e_starts, e_ends = [], []
    for e_st, e_en in exList:
        e_starts.append(e_st)
        e_ends.append(e_en)

    startAll, endAll = min(e_starts), max(e_ends)

    gene_rows = ann_df.loc[(ann_df["gene"].str.upper() == geneID.strip().upper()) & (ann_df["chrom"] == chrom)]
    exonStartList = {}
    exonEndList = {}
    exonCountList = {}

    isoforms = 0
    for index, row in gene_rows.iterrows():
        exonCount = row["exonCount"]
        exonCountList[isoforms] = exonCount
        exonStarts = list(filter(None, row["exonStarts"].split(",")))
        exonEnds = list(filter(None, row["exonEnds"].split(",")))
        exonStartList[isoforms] = exonStarts
        exonEndList[isoforms] = exonEnds
        isoforms += 1

    title = "" + chrom + ":" + start + "-" + end + "(" + geneID + ")"
    print("Number of isoforms found for gene " + geneID + ": " + str(isoforms))
    number_of_subplots = len(s1_namelist) + len(s2_namelist) + 1
    fig, axes = plt.subplots(nrows=number_of_subplots, ncols=1)

    print("Generating read coverage plots ...")
    for i, ax1 in enumerate(axes[0 : len(s1_namelist)]):
        y_limit = Generate_read_coverate_plot(
            ax1,
            input1_dir,
            s1_namelist[i],
            chrom,
            geneID,
            int(start),
            int(end),
            startAll,
            endAll,
            g1_name,
            1,
        )
        if i == 0:
            ax1.set_title(title, color="black", fontsize=14)
    for i, ax2 in enumerate(axes[len(s1_namelist) : number_of_subplots - 1]):
        j = i - len(s1_namelist)
        y_limit = Generate_read_coverate_plot(
            ax2,
            input2_dir,
            s2_namelist[j],
            chrom,
            geneID,
            int(start),
            int(end),
            startAll,
            endAll,
            g2_name,
            2,
        )

    print("Generating annotation plots ...")
    ax3 = axes[number_of_subplots - 1]
    Generate_annotation_plot(
        ax3,
        isoforms,
        exonCountList,
        exonStartList,
        exonEndList,
        int(start),
        int(end),
        startAll,
        endAll,
    )
    ax3.set_ylabel("Annotation", fontsize=12)
    ax3.set_xlabel("Position", fontsize=12)
    ax3.spines["top"].set_color("none")
    ax3.spines["left"].set_color("none")
    ax3.spines["right"].set_color("none")
    ax3.tick_params(labelcolor="w", top=False, left=False, right=False)

    plt.savefig(os.path.join(output_dir, title + ".png"))
    plt.savefig(os.path.join(output_dir, title + ".eps"), format="eps", dpi=100)
    print(title + " plotted successfully.")
    y_limit = 0


######### Main program ###########
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

# Read in the configuration values
config = configparser.ConfigParser()
config.read("./config/config.ini")

# Input values
gtf_file = config["INPUT_FILES"]["gtf_file"]
ref_file = config["INPUT_FILES"]["ref_file"]
group1_folder = config["INPUT_FILES"]["group1_folder"]
group2_folder = config["INPUT_FILES"]["group2_folder"]
annotation_file = config["INPUT_FILES"]["annotation_file"]

# Output values
as_events_output_dir = config["OUTPUT"]["as_events_output_dir"]
as_quant_output_dir = config["OUTPUT"]["as_quant_output_dir"]
gtf_output_file = config["OUTPUT"]["gtf_output_file"]

# Parameters
samtools = config["PARAMETERS"]["samtools_dir"]
cores = config["PARAMETERS"]["cores"]
species = config["PARAMETERS"]["species"]
method = config["PARAMETERS"]["method"]
novel = config["PARAMETERS"]["novel"]


output_dir = as_quant_output_dir
input1_dir = group1_folder
input2_dir = group2_folder

g1_name, g2_name = os.path.basename(input1_dir), os.path.basename(input2_dir)

if species == "hg38" or species == "hg19":
    chromosomes = chromosomes_h
elif species == "mm10":
    chromosomes = chromosomes_m

s1_namelist, s2_namelist = list_dirs(input1_dir), list_dirs(input2_dir)

ann_df = pd.read_csv(annotation_file, delimiter="\t")
print("Collecting annotation and input data for AS-Quant run ...")
#### convert the whole annotation into a dictionary for faster use
ChromDict = MakeFullDictionary(ann_df, chromosomes)
#### merge the exons intervals #####
ChromDict_merged = merge_ChromDict(ChromDict, chromosomes)

while True:
    region = input("Enter the range: (chr:gene:start-end): ")
    if region.lower() == "exit":
        print("Exiting from AS-Quant plots generation.")
        break
    Process_user_inputs(
        region,
        ChromDict_merged,
        input1_dir,
        input2_dir,
        s1_namelist,
        s2_namelist,
        g1_name,
        g2_name,
        output_dir,
        ann_df,
    )
    # 👇 Exit when user put 'exit' and presses Enter
