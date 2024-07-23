import os
import time
from methods.print_and_log import print_and_log


# This function determines if items are unique in a list
def unique(input_list):
    uniques = []
    for item in input_list:
        if item not in uniques:
            uniques.append(item)
    return uniques


# This function finds alternative splicing events in the GTF file
def find_as_events(gtf, prefix, output_dir, logger):
    # Step 1: Find alternative splicing events
    print_and_log("----------------------------------------------------------------", logger)
    print_and_log("| STEP 1: Finding alternative splicing events...               |", logger)
    print_and_log("----------------------------------------------------------------", logger)

    # Start logging time and create the output directories
    start = time.perf_counter()
    gtf_input_file = open(gtf)
    output_dir = os.path.join(output_dir, prefix)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Create the output files for the alternative splicing events
    a3ss_output = open(os.path.join(output_dir, "A3SS.csv"), "w")
    a5ss_output = open(os.path.join(output_dir, "A5SS.csv"), "w")
    se_output = open(os.path.join(output_dir, "SE.csv"), "w")
    mxe_output = open(os.path.join(output_dir, "MXE.csv"), "w")
    afe_output = open(os.path.join(output_dir, "AFE.csv"), "w")
    ale_output = open(os.path.join(output_dir, "ALE.csv"), "w")
    ri_output = open(os.path.join(output_dir, "RI.csv"), "w")

    # Create the initial headers for the output files
    ce_header = "chrom\tgene\texonStart\texonEnd\tstrand"
    mxe_header = "chrom\tgene\texon1Start\texon1End\texon2Start\texon2End\tstrand"
    alt_ss_header = "chrom\tgene\tlongExonStart\tlongExonEnd\tshortExonStart\tshortExonEnd\tstrand"
    alt_fl_headerder = "chrom\tgene\tdistalExonStart\tdistalExonEnd\tproximalExonStart\tproximalExonEnd\tstrand"
    ri_header = "chrom\tgene\texonStart\tExonEnd\tstrand"

    # Write headers to the output files
    se_output.write(ce_header + "\n")
    mxe_output.write(mxe_header + "\n")
    a3ss_output.write(alt_ss_header + "\n")
    a5ss_output.write(alt_ss_header + "\n")
    afe_output.write(alt_fl_headerder + "\n")
    ale_output.write(alt_fl_headerder + "\n")
    ri_output.write(ri_header + "\n")

    # Process the GTF file
    chunk_size = 1000
    gene_group = {}
    genes = {}
    supplement = {}

    # Iterate through the GTF file
    for line in gtf_input_file:
        # Extract the chrom, transcript_type, start_pos, end_pos, group, and strand from the line
        elements = line.strip().split("\t")
        chrom = elements[0]
        transcript_type = elements[2]
        start_pos = elements[3]
        end_pos = elements[4]
        group = range(int(int(start_pos) / chunk_size), int(int(end_pos) / chunk_size) + 1)  ## groups this line could belong to
        group = list(set(group))
        strand = elements[6]

        # Desc holds the gene_id, transcript_id, and gene_name
        desc = elements[8].split(";")
        gene_id = ["", ""]
        transcript_id = ["", ""]
        gene_name = ["", "NA"]

        # Extract and update the gene_id, transcript_id, and gene_name
        for desc_element in desc:
            # Guard for the last element
            if len(desc_element.strip()) < 2 or len(desc_element.strip().split(" ")) < 2:
                continue

            # Extract the key value pair from the description element
            desc_key = desc_element.strip().split(" ")[0]
            desc_val = desc_element.strip().split(" ")[1]
            if desc_key.upper() == "GENE_ID":
                gene_id = [desc_key, desc_val]
            elif desc_key.upper() == "TRANSCRIPT_ID":
                transcript_id = [desc_key, desc_val]
            elif desc_key.upper() == "GENE_NAME":
                gene_name = [desc_key, desc_val]

        # Check if the gene_id and transcript_id are correct
        if gene_id[0].upper() != "GENE_ID" or transcript_id[0].upper() != "TRANSCRIPT_ID":
            logger.debug("This line does not have correct description for gene_id or transcript_id: %s, %s" % (gene_id, transcript_id))
            logger.debug("Incorrect description: %s" % elements)
            continue

        # Update the gene_group, genes, and supplement dictionaries
        for i in group:
            if i in gene_group:
                gene_group[i].append(gene_id[1])
            else:
                gene_group[i] = [gene_id[1]]

        # Add exon to the gene and appropriate transcript
        if transcript_type == "exon":
            # Check if we've already processed the gene and transcript
            if gene_id[1] in genes:
                if transcript_id[1] in genes[gene_id[1]]:
                    # Append exon to the transcript
                    genes[gene_id[1]][transcript_id[1]].append([int(start_pos), int(end_pos)])
                else:
                    # Add the first exon to the transcript
                    genes[gene_id[1]][transcript_id[1]] = [[int(start_pos), int(end_pos)]]  ## add first exon
            else:
                # Add the first exon to the gene
                genes[gene_id[1]] = {}
                genes[gene_id[1]][transcript_id[1]] = [[int(start_pos), int(end_pos)]]  ## add first exon
                supplement[gene_id[1]] = [gene_name[1], chrom, strand]  ## geneName, chromosom and strand

    # Remove duplicates from the gene_group
    for group in gene_group:
        gene_group[group] = list(set(gene_group[group]))

    # Log the gene statistics before procesing
    logger.info("~~~~~~~~~~~~~~~~~~~~~~~~~~ GENE STATS ~~~~~~~~~~~~~~~~~~~~~~~~~~")
    num_genes = len(genes)
    num_transcripts = 0
    num_one_transcripts = 0
    num_exons = 0
    num_one_exon = 0
    num_one_transcript_exon = 0

    # Gather gene statistics
    for gene in genes:  ##
        # Update the gene count
        num_transcripts += len(genes[gene])

        # Check if the gene has only one transcript
        if len(genes[gene]) == 1:
            num_one_transcripts += 1

        # Iterate through the transcripts
        for transcript in genes[gene]:
            # Update the exon count
            num_exons += len(genes[gene][transcript])

            # Check if the transcript has only one exon
            if len(genes[gene][transcript]) == 1:
                num_one_exon += 1

                # Check if the gene has only one transcript with one exon
                if len(genes[gene]) == 1:
                    num_one_transcript_exon += 1

    # Log the gene statistics
    logger.info("%d distinct gene ID in the gtf file" % num_genes)
    logger.info("%d distinct transcript ID in the gtf file" % num_transcripts)
    logger.info("%d one-transcript genes in the gtf file" % num_one_transcripts)
    logger.info("%d exons in the gtf file" % num_exons)
    logger.info("%d one-exon transcripts in the gtf file" % num_one_exon)
    logger.info("%d one-transcript genes with only one exon in the transcript" % num_one_transcript_exon)

    # Avoid division by zero
    if num_genes > 0:
        logger.info("Average number of transcripts per gene is %f" % (float(num_transcripts) / num_genes))

    # Avoid division by zero
    if num_transcripts > 0:
        logger.info("Average number of exons per transcript is %f" % (float(num_exons) / num_transcripts))

    # Avoid division by zero
    if (num_transcripts - num_one_exon) > 0:  ## to avoid divided by zero exception
        logger.info("Average number of exons per transcript excluding one-exon tx is %f" % (float(num_exons - num_one_exon) / (num_transcripts - num_one_exon)))

    logger.info("~~~~~~~~~~~~~~~~~~~~~~~ GENE GROUP STATS ~~~~~~~~~~~~~~~~~~~~~~~")
    tgi = 0
    for group in gene_group:
        tgi += len(gene_group[group])
    logger.info("There are total of %d groups and %d genes in gene_group" % (len(gene_group), tgi))
    logger.info("The average number of genes in each group is %f" % (float(tgi) / len(gene_group)))
    logger.info("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

    # Sort the gene transcripts
    for gene_id in genes:
        for transcript_id in genes[gene_id]:
            if len(genes[gene_id][transcript_id]) == 1:
                continue
            genes[gene_id][transcript_id] = sorted(genes[gene_id][transcript_id])

    # This function determines if an exon is fully contained in an internal exon
    def fullyContainedInInternalExon(myExon, myGeneID):
        # Iterate through the transcripts
        for mytranscript_id in genes[myGeneID]:
            # For each internal exon
            for intExon in genes[myGeneID][mytranscript_id][1:-1]:
                # Check if the exon is fully contained
                if intExon[0] <= myExon[0] and intExon[1] >= myExon[1]:  ## fully contained
                    return True
        return False

    # Dictionaries for alternative splicing events
    s_events = {}
    mx_events = {}
    ss3_events = {}
    ss5_events = {}
    afe_eventss = {}
    ale_events = {}
    ri_events = {}

    # Accumulator values
    num_skipping_events = 0
    num_mxe_events = 0
    num3 = 0
    num5 = 0
    num_afe = 0
    num_ale = 0
    num_ri = 0

    # Accumulator values for duplicates
    duplicate_se = 0
    duplicate_mxe = 0
    duplicate_ss3 = 0
    duplicate_ss5 = 0
    duplicate_ri = 0

    # Iterate through each of the gene IDs
    for gene_id in genes:  ## process each gene
        supInfo = supplement[gene_id]  ## supplementmentary info
        if len(genes[gene_id]) == 1:  ## only one transcript, alt SS event is imposible
            continue  ## go to the next geneID
        else:  ## look for alt SS event
            de = {}  ## distinct exons
            for tID in genes[gene_id]:  ## sort each transcript, merge and get distinct exons
                if len(genes[gene_id][tID]) == 1:  ## only one exon, skip it
                    continue  ## next transcript..
                genes[gene_id][tID] = sorted(genes[gene_id][tID])  ## sort each transcript
                for exon in genes[gene_id][tID]:
                    de[exon[0], exon[1]] = 1  ## it's okay to overwrite

            for ce in de:  ## for each exon in distinct exon dictionary
                uf = []  ## upstream flanking exons
                df = []  ## downstream flanking exons
                for tID in genes[gene_id]:  ## examine each transcript to see if it contains the given exon
                    if [ce[0], ce[1]] in genes[gene_id][tID]:  ## this exon is in the transcript
                        eInd = genes[gene_id][tID].index([ce[0], ce[1]])
                        if 0 < eInd:  ## it definitely has upstream flanking exon
                            uf.append(genes[gene_id][tID][eInd - 1])
                        if eInd < len(genes[gene_id][tID]) - 1:  ## it definitely has downstream flanking exon
                            df.append(genes[gene_id][tID][eInd + 1])

                ## getting uniq upstream flanking exons and downstream flanking exons
                uf = unique(uf)
                df = unique(df)
                #### exon skipping (cassette exon) events ###
                for i in range(0, len(uf)):  ### going through the upstream flanking exons
                    f1 = uf[i]  ## first flanking exon
                    for j in range(0, len(df)):  ## going through the downstream falnking exons to see if there is skipping junction
                        f2 = df[j]
                        for tID in genes[gene_id]:  ## examine each transcript, to see if it has an edge connecting f1 and f2
                            if len(genes[gene_id][tID]) < 2:  ## less than two exons, skip it
                                continue

                            ########################## not requiring exact falnking exons ###########
                            for i in range(0, len(genes[gene_id][tID]) - 1):  ## for each exon in the tx
                                e_1 = genes[gene_id][tID][i]
                                e_2 = genes[gene_id][tID][i + 1]
                                if e_1[1] == f1[1] and e_2[0] == f2[0]:  ## this tx has an edge connecting f1 and f2 but does not have ce
                                    key = supInfo[1] + ":" + str(ce[0] - 1) + ":" + str(ce[1]) + ":" + str(f1[1]) + ":" + str(f2[0] - 1)  ## target exon and skipping junction
                                    if key in s_events:  ## already have this skipping events
                                        duplicate_se += 1
                                        continue  ## next transcript
                                    else:  ## new key, write it
                                        s_events[key] = 1
                                        num_skipping_events += 1
                                        se_output.write(supInfo[1] + "\t" + gene_id + "\t" + str(f1[0] - 1) + "\t" + str(f1[1]) + "\t" + supInfo[2] + "\n")

                ##### mutually exclusive events #####
                for i in range(0, len(uf)):  ### going through the upstream flanking exons
                    f1 = uf[i]  ## first flanking exon
                    for j in range(0, len(df)):  ## going through the downstream falnking exons to see if there is skipping junction
                        f2 = df[j]
                        for tID in genes[gene_id]:  ## examine each transcript
                            if [ce[0], ce[1]] in genes[gene_id][tID] or f1 not in genes[gene_id][tID] or f2 not in genes[gene_id][tID]:  ## ce in, f1 or f2 is not in..do not examine
                                continue  ### go to next transcript in genes[gene_id]
                            else:  ## this transcript does not have ce and has both f1 and f2, let's take a look
                                fromF1 = genes[gene_id][tID].index(f1)
                                fromF2 = genes[gene_id][tID].index(f2)
                                mxe = genes[gene_id][tID][fromF1 + 1]  ## candidate mxe
                                if (fromF1 + 1 == fromF2 - 1) and (mxe[0] > ce[1]):  ### this exon is the right one
                                    goodMXE = True
                                    if goodMXE:  ### it's okay to write out
                                        # key = supInfo[1]+':'+str(ce[0]-1)+':'+str(ce[1])+':'+str(mxe[0]-1)+':'+str(mxe[1])+':'+str(f1[0]-1)+':'+str(f1[1])+':'+str(f2[0]-1)+':'+str(f2[1]);
                                        key = supInfo[1] + ":" + str(ce[0] - 1) + ":" + str(ce[1]) + ":" + str(mxe[0] - 1) + ":" + str(mxe[1]) + ":" + str(f1[1]) + ":" + str(f2[0] - 1)
                                        if key in mx_events:  ## duplicate,
                                            duplicate_mxe += 1
                                        else:
                                            num_mxe_events += 1
                                            mx_events[key] = 1
                                            mxe_output.write(supInfo[1] + "\t" + gene_id + "\t" + str(mxe[0] - 1) + "\t" + str(mxe[1]) + "\t" + str(f1[0] - 1) + "\t" + str(f1[1]) + "\t" + supInfo[2] + "\n")
                                            # if tID.find("novel") >= 0:  # It involves novel junction
                                            #     neFile_mxe.write(str(num_mxe_events) + "\t" + gene_id + "\t" + "\t".join(supInfo) + "\t" + str(ce[0] - 1) + "\t" + str(ce[1]) + "\t" + str(mxe[0] - 1) + "\t" + str(mxe[1]) + "\t" + str(f1[0] - 1) + "\t" + str(f1[1]) + "\t" + str(f2[0] - 1) + "\t" + str(f2[1]) + "\n")

                #### alt-3 and alt-5 events ###
                for i in range(0, len(uf) - 1):  ### going through the upstream flanking exons
                    e = uf[i]
                    for j in range(i + 1, len(uf)):
                        u = uf[j]
                        if e[0] == u[0]:  ## it is alt SS event, because uf is derived from SET
                            # key = supInfo[1]+':'+str(ce[0]-1)+':'+str(ce[1])+':'+str(e[0]-1)+':'+str(max(e[1],u[1]))+':'+str(e[0]-1)+':'+str(min(e[1],u[1]));
                            key = supInfo[1] + ":" + str(min(e[1], u[1])) + ":" + str(max(e[1], u[1])) + ":" + str(ce[0] - 1)
                            if supInfo[2] == "+":  ## positive strand. alt-5 event
                                if key in ss5_events:  ## duplicate
                                    duplicate_ss5 += 1
                                else:
                                    ss5_events[key] = 1
                                    num5 += 1
                                    a5ss_output.write(supInfo[1] + "\t" + gene_id + "\t" + str(e[0] - 1) + "\t" + str(min(e[1], u[1])) + "\t" + str(ce[0] - 1) + "\t" + str(ce[1]) + "\t" + supInfo[2] + "\n")
                            else:  ## neg strand. alt-3 event
                                if key in ss3_events:  ## duplicate
                                    duplicate_ss3 += 1
                                else:
                                    ss3_events[key] = 1
                                    num3 += 1
                                    a3ss_output.write(supInfo[1] + "\t" + gene_id + "\t" + str(e[0] - 1) + "\t" + str(min(e[1], u[1])) + "\t" + str(ce[0] - 1) + "\t" + str(ce[1]) + "\t" + supInfo[2] + "\n")

                for i in range(0, len(df) - 1):  ### going through the downstream flanking exons
                    e = df[i]
                    for j in range(i + 1, len(df)):
                        d = df[j]
                        if e[1] == d[1]:  ## it is alt SS event, because uf is derived from SET
                            key = supInfo[1] + ":" + str(ce[0] - 1) + ":" + str(ce[1]) + ":" + str(min(e[0], d[0]) - 1) + ":" + str(e[1]) + ":" + str(max(e[0], d[0]) - 1) + ":" + str(e[1])
                            key = supInfo[1] + ":" + str(ce[1]) + ":" + str(min(e[0], d[0]) - 1) + ":" + str(max(e[0], d[0]) - 1)
                            if supInfo[2] == "+":  ## positive strand. alt-3 event
                                if key in ss3_events:  ## duplicate
                                    duplicate_ss3 += 1
                                else:
                                    ss3_events[key] = 1
                                    num3 += 1
                                    a3ss_output.write(supInfo[1] + "\t" + gene_id + "\t" + str(max(e[0], d[0]) - 1) + "\t" + str(e[1]) + "\t" + str(ce[0] - 1) + "\t" + str(ce[1]) + "\t" + supInfo[2] + "\n")
                                    # if tID.find("novel") >= 0:  # It involves novel junction
                                    #     neFile_3.write(str(num3) + "\t" + gene_id + "\t" + "\t".join(supInfo) + "\t" + str(min(e[0], d[0]) - 1) + "\t" + str(e[1]) + "\t" + str(max(e[0], d[0]) - 1) + "\t" + str(e[1]) + "\t" + str(ce[0] - 1) + "\t" + str(ce[1]) + "\n")
                            else:  ## neg strand. alt-5 event
                                if key in ss5_events:  ## duplicate
                                    duplicate_ss5 += 1
                                else:
                                    ss5_events[key] = 1
                                    num5 += 1
                                    a5ss_output.write(supInfo[1] + "\t" + gene_id + "\t" + str(max(e[0], d[0]) - 1) + "\t" + str(e[1]) + "\t" + str(ce[0] - 1) + "\t" + str(ce[1]) + "\t" + supInfo[2] + "\n")
                                    # if tID.find("novel") >= 0:  # It involves novel junction
                                    #     neFile_5.write(str(num5) + "\t" + gene_id + "\t" + "\t".join(supInfo) + "\t" + str(min(e[0], d[0]) - 1) + "\t" + str(e[1]) + "\t" + str(max(e[0], d[0]) - 1) + "\t" + str(e[1]) + "\t" + str(ce[0] - 1) + "\t" + str(ce[1]) + "\n")

                ### Alternative First Exon and Alternative Last Exon

                for tID in genes[gene_id]:  ## examine each transcript to see if the given exon is in a given transcript (index should be 1 or len-2)
                    tLen = len(genes[gene_id][tID])  ## length of a transcript
                    if tLen < 2:  ## not enough exons, skip this
                        continue  ## process next transcript

                    if [ce[0], ce[1]] == genes[gene_id][tID][1]:  ## current exon is the 2nd exon in the tx
                        fEx = genes[gene_id][tID][0]  ## firstExon, find other tx with different fEX (non overlapping)
                        if fullyContainedInInternalExon(fEx, gene_id):  ## fEx is fully contained.. process next transcript
                            continue
                        for ctranscript_id in genes[gene_id]:  ## need to examine each tx, candidate transcript id
                            if len(genes[gene_id][ctranscript_id]) < 2:  ## not enough exons, skip this
                                continue  ## process next candidate transcript
                            if [ce[0], ce[1]] == genes[gene_id][ctranscript_id][1]:  ## the target exon is the 2nd exon in the ctx
                                cfEx = genes[gene_id][ctranscript_id][0]  ## candidate first exon
                                if cfEx[0] > fEx[1] or fEx[0] > cfEx[1]:  ### non-overlapping exon with smaller coord
                                    #### should not fully contained in an internal exon of other transcript
                                    if fullyContainedInInternalExon(cfEx, gene_id):  ## cfEx is fully contained.. process next candidate transcript
                                        continue

                                    if supInfo[2] == "+":  ## positive strand. AFE, alt first exon event
                                        key = supInfo[1] + ":" + str(ce[0] - 1) + ":" + str(ce[1]) + ":"
                                        key += str(min(fEx[0], cfEx[0]) - 1) + ":" + str(min(fEx[1], cfEx[1])) + ":"
                                        key += str(max(fEx[0], cfEx[0]) - 1) + ":" + str(max(fEx[1], cfEx[1]))
                                        if key in afe_eventss:  ## already have this one..
                                            pass  ## do nothing
                                        else:  ## new AFE
                                            afe_eventss[key] = 1
                                            num_afe += 1
                                            afe_output.write(supInfo[1] + "\t" + gene_id + "\t" + str(max(fEx[0], cfEx[0]) - 1) + "\t" + str(max(fEx[1], cfEx[1])) + "\t" + str(ce[0] - 1) + "\t" + str(ce[1]) + "\t" + supInfo[2] + "\n")
                                    else:  ## neg strand. ALE, alt last exon event
                                        key = supInfo[1] + ":" + str(ce[0] - 1) + ":" + str(ce[1]) + ":"
                                        key += str(min(fEx[0], cfEx[0]) - 1) + ":" + str(min(fEx[1], cfEx[1])) + ":"
                                        key += str(max(fEx[0], cfEx[0]) - 1) + ":" + str(max(fEx[1], cfEx[1]))
                                        if key in ale_events:  ## already have this one..
                                            pass  ## do nothing
                                        else:  ## new ALE
                                            ale_events[key] = 1
                                            num_ale += 1
                                            ale_output.write(supInfo[1] + "\t" + gene_id + "\t" + str(max(fEx[0], cfEx[0]) - 1) + "\t" + str(max(fEx[1], cfEx[1])) + "\t" + str(ce[0] - 1) + "\t" + str(ce[1]) + "\t" + supInfo[2] + "\n")

                    if [ce[0], ce[1]] == genes[gene_id][tID][-2]:  ## current exon is the 2nd to the last exon in the tx
                        fEx = genes[gene_id][tID][-1]  ## lastExon, find other tx with different fEX (non overlapping)
                        if fullyContainedInInternalExon(fEx, gene_id):  ## fEx is fully contained.. process next transcript
                            continue
                        for ctranscript_id in genes[gene_id]:  ## need to examine each tx, candidate transcript id
                            if len(genes[gene_id][ctranscript_id]) < 2:  ## not enough exons, skip this
                                continue  ## process next candidate transcript
                            if [ce[0], ce[1]] == genes[gene_id][ctranscript_id][-2]:  ## the target exon is the 2nd exon in the ctx
                                cfEx = genes[gene_id][ctranscript_id][-1]  ## candidate last exon
                                if cfEx[0] > fEx[1] or fEx[0] > cfEx[1]:  ### non-overlapping exon with smaller coord
                                    #### should not fully contained in an internal exon of other transcript
                                    if fullyContainedInInternalExon(cfEx, gene_id):  ## cfEx is fully contained.. process next candidate transcript
                                        continue

                                    if supInfo[2] == "-":  ## negative strand. AFE, alt first exon event
                                        key = supInfo[1] + ":" + str(ce[0] - 1) + ":" + str(ce[1]) + ":"
                                        key += str(min(fEx[0], cfEx[0]) - 1) + ":" + str(min(fEx[1], cfEx[1])) + ":"
                                        key += str(max(fEx[0], cfEx[0]) - 1) + ":" + str(max(fEx[1], cfEx[1]))
                                        if key in afe_eventss:  ## already have this one..
                                            pass  ## do nothing
                                        else:  ## new AFE
                                            afe_eventss[key] = 1
                                            num_afe += 1
                                            afe_output.write(supInfo[1] + "\t" + gene_id + "\t" + str(min(fEx[0], cfEx[0]) - 1) + "\t" + str(min(fEx[1], cfEx[1])) + "\t" + str(ce[0] - 1) + "\t" + str(ce[1]) + "\t" + supInfo[2] + "\n")
                                    else:  ## pos strand. ALE, alt last exon event
                                        key = supInfo[1] + ":" + str(ce[0] - 1) + ":" + str(ce[1]) + ":"
                                        key += str(min(fEx[0], cfEx[0]) - 1) + ":" + str(min(fEx[1], cfEx[1])) + ":"
                                        key += str(max(fEx[0], cfEx[0]) - 1) + ":" + str(max(fEx[1], cfEx[1]))
                                        if key in ale_events:  ## already have this one..
                                            pass  ## do nothing
                                        else:  ## new ALE
                                            ale_events[key] = 1
                                            num_ale += 1
                                            ale_output.write(supInfo[1] + "\t" + gene_id + "\t" + str(min(fEx[0], cfEx[0]) - 1) + "\t" + str(min(fEx[1], cfEx[1])) + "\t" + str(ce[0] - 1) + "\t" + str(ce[1]) + "\t" + supInfo[2] + "\n")

                ### Retained Intron events
                for i in range(0, len(uf)):  ### going through the upstream flanking exons
                    f1 = uf[i]  ## first flanking exon
                    for tID in genes[gene_id]:  ## examine each transcript
                        if [f1[0], ce[1]] in genes[gene_id][tID]:  ## there is an exon starts from f1 ends at ce, it is retained intron
                            # key=supInfo[1]+':'+str(f1[0]-1)+':'+str(ce[1])+':'+str(f1[0]-1)+':'+str(f1[1])+':'+str(ce[0]-1)+':'+str(ce[1]);
                            key = supInfo[1] + ":" + str(f1[1]) + ":" + str(ce[0] - 1)
                            if key in ri_events:  ## already have this skipping events
                                duplicate_ri += 1
                                continue  ## next transcript
                            else:  ## new key, write it
                                ri_events[key] = 1
                                num_ri += 1
                                ri_output.write(supInfo[1] + "\t" + gene_id + "\t" + str(f1[0] - 1) + "\t" + str(f1[1]) + "\t" + supInfo[2] + "\n")
                                # if tID.find("novel") >= 0:  # It involves novel junction
                                #     neFile_ri.write(str(num_ri) + "\t" + gene_id + "\t" + "\t".join(supInfo) + "\t" + str(f1[0] - 1) + "\t" + str(ce[1]) + "\t" + str(f1[0] - 1) + "\t" + str(f1[1]) + "\t" + str(ce[0] - 1) + "\t" + str(ce[1]) + "\n")

                for i in range(0, len(df)):  ### going through the downstream flanking exons
                    f1 = df[i]  ## first flanking exon
                    for tID in genes[gene_id]:  ## examine each transcript
                        if [ce[0], f1[1]] in genes[gene_id][tID]:  ## there is an exon starts from ce ends at f1, it is retained intron
                            key = supInfo[1] + ":" + str(ce[1]) + ":" + str(f1[0] - 1)
                            if key in ri_events:  ## already have this skipping events
                                duplicate_ri += 1
                                continue  ## next transcript
                            else:  ## new key, write it
                                ri_events[key] = 1
                                num_ri += 1
                                ri_output.write(supInfo[1] + "\t" + gene_id + "\t" + str(f1[0] - 1) + "\t" + str(f1[1]) + "\t" + supInfo[2] + "\n")

    # Close all the output files
    gtf_input_file.close()
    a3ss_output.close()
    a5ss_output.close()
    se_output.close()
    mxe_output.close()
    afe_output.close()
    ale_output.close()
    ri_output.close()

    # Log the completion status and statistics
    total_events = num_skipping_events + num_mxe_events + num5 + num3 + num_afe + num_ale + num_ri
    print_and_log(f"{num_skipping_events} Exon skipping events ({num_skipping_events/total_events*100:.2f}%)", logger)
    print_and_log(f"{num_mxe_events} MX events ({num_mxe_events/total_events*100:.2f}%)", logger)
    print_and_log(f"{num5} Alt 5 SS events ({num5/total_events*100:.2f}%)", logger)
    print_and_log(f"{num3} Alt 3 SS events ({num3/total_events*100:.2f}%)", logger)
    print_and_log(f"{num_afe} AFE events ({num_afe/total_events*100:.2f}%)", logger)
    print_and_log(f"{num_ale} ALE events ({num_ale/total_events*100:.2f}%)", logger)
    print_and_log(f"{num_ri} RI events ({num_ri/total_events*100:.2f}%)", logger)
    logger.info(f"{duplicate_se} duplicate skipping events")
    logger.info(f"{duplicate_mxe} duplicate MXE events")
    logger.info(f"{duplicate_ss5} duplicate alt-5 SS events")
    logger.info(f"{duplicate_ss3} duplicate alt-3 SS events")
    logger.info(f"{duplicate_ri} duplicate RI events")
    print_and_log(f"Step 1 Completed in {time.perf_counter() - start:.2f} seconds", logger)
    print_and_log("----------------------------------------------------------------\n", logger)
    logger.handlers[0].stream.write("\n\n")
