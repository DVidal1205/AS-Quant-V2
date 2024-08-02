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
    se_header = "chrom\tgene\texonStart\texonEnd\tstrand"
    mxe_header = "chrom\tgene\texon1Start\texon1End\texon2Start\texon2End\tstrand"
    alt_ss_header = "chrom\tgene\tlongExonStart\tlongExonEnd\tshortExonStart\tshortExonEnd\tstrand"
    alt_fl_headerder = "chrom\tgene\tdistalExonStart\tdistalExonEnd\tproximalExonStart\tproximalExonEnd\tstrand"
    ri_header = "chrom\tgene\texonStart\texonEnd\tstrand"

    # Write headers to the output files
    se_output.write(se_header + "\n")
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

    # Data stucture representation:
    #
    #
    # gene_group:
    # {
    #     0: ["gene1", "gene2"],
    #     1: ["gene2", "gene3"],
    #     2: ["gene3"]
    # }
    #
    # genes:
    # {
    #     "gene1": {
    #         "transcript1": [[100, 200], [300, 400]],
    #         "transcript2": [[150, 250], [350, 450]]
    #     },
    #     "gene2": {
    #         "transcript3": [[500, 600]],
    #         "transcript4": [[700, 800], [900, 1000]]
    #     },
    #     "gene3": {
    #         "transcript5": [[1100, 1200]]
    #     }
    # }
    #
    # supplement:
    # {
    #     "gene1": ["geneName1", "chr1", "+"],
    #     "gene2": ["geneName2", "chr2", "-"],
    #     "gene3": ["geneName3", "chr3", "+"]
    # }

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

    # Dictionaries for alternative splicing events
    s_events = {}
    mx_events = {}
    ss3_events = {}
    ss5_events = {}
    afe_events = {}
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

    # Iterate through each of the gene IDs, then each of the transcripts, and finally each of the exons
    for gene_id in genes:
        sup_info = supplement[gene_id]

        # Check if the gene has only one transcript, meaning no alternative splicing events are possible
        if len(genes[gene_id]) == 1:
            continue

        # DISTINCT EXONS DICTIONARY
        # key: start and end of the exon, as a tuple
        # value: 1, to indicate that the exon is distinct
        distinct_exons = {}

        # Iterate through each of the transcripts in the gene, populating the distinct exons dictionary
        for transcript_id in genes[gene_id]:
            # Skip if there is only one exon
            if len(genes[gene_id][transcript_id]) == 1:  ## only one exon, skip it
                continue

            # Sort each of the transcripts in the gene
            genes[gene_id][transcript_id] = sorted(genes[gene_id][transcript_id])

            # Iterate through each of the exons in the gene
            for exon in genes[gene_id][transcript_id]:
                # Add the exon to the distinct exons dictionary by using the start and end values as keys
                distinct_exons[exon[0], exon[1]] = 1

        # Iterate through each of the distinct exons
        for current_exon in distinct_exons:
            # Initialize the upstream and downstream lists
            upstream = []
            downstream = []

            # Iterate through each of the transcripts in the current gene
            for transcript_id in genes[gene_id]:
                # Check if the current exon is in the current transcript
                if [current_exon[0], current_exon[1]] in genes[gene_id][transcript_id]:
                    # Extract the index of the current exon in the transcript
                    exon_index = genes[gene_id][transcript_id].index([current_exon[0], current_exon[1]])

                    # If the exon index is greater than 0, then it must have an upstream flanking exon
                    if exon_index > 0:
                        upstream.append(genes[gene_id][transcript_id][exon_index - 1])

                    # If the exon index is less than the length of the gene - 1, then it must have a downstream flanking exon
                    if exon_index < len(genes[gene_id][transcript_id]) - 1:
                        downstream.append(genes[gene_id][transcript_id][exon_index + 1])

            # Only consider unique exons for the upstream and downstream lists
            upstream = unique(upstream)
            downstream = unique(downstream)

            """
            1. SKIPPING EXON EVENTS
            
            This section of code processes exon skipping (SE) events for a given gene. The goal is to identify cases where an exon 
            is skipped in some transcripts of the gene. It does this by examining the relationships between upstream and downstream exons 
            flanking the exon of interest.
            
            This code iterates through all possible pairs of upstream and downstream exons, checking if the exon of interest is present in
            the current transcript. If there is an edge connecting the first and second flanking exons but not the exon of interest, then a
            skipping event is detected.
            """

            # Iterate through the upstream exons
            for i in range(0, len(upstream)):
                # Extract the first flanking exon
                first_flanking = upstream[i]
                # Iterate through the downstream exons
                for j in range(0, len(downstream)):
                    # Extract the second flanking exon
                    second_flanking = downstream[j]

                    # Iterate through the transcripts in the gene
                    for transcript_id in genes[gene_id]:
                        # If there are less than two exons, skip it
                        if len(genes[gene_id][transcript_id]) < 2:
                            continue

                        # Iterate through neighboring exon pairs
                        for i in range(0, len(genes[gene_id][transcript_id]) - 1):
                            # Extract the current exon and the next exon
                            exon_1 = genes[gene_id][transcript_id][i]
                            exon_2 = genes[gene_id][transcript_id][i + 1]

                            # If this transcript has an edge connecting the first and second flanking exons but does not have the current exon
                            if exon_1[1] == first_flanking[1] and exon_2[0] == second_flanking[0]:
                                # Create a unique key for the skipping event
                                key = f"{sup_info[1]}:{current_exon[0] - 1}:{current_exon[1]}:{first_flanking[1]}:{second_flanking[0] - 1}"

                                # If the key is already in the skipping events dictionary, increment the duplicate count and continue
                                if key in s_events:
                                    duplicate_se += 1
                                    continue

                                # Otherwise, add the key to the skipping events dictionary and increment the skipping event count
                                else:
                                    s_events[key] = 1
                                    num_skipping_events += 1
                                    se_output.write(f"{sup_info[1]}\t{gene_id}\t{first_flanking[0] - 1}\t{first_flanking[1]}\t{sup_info[2]}\n")

            """
            2. MUTUALLY EXCLUSIVE EXON EVENTS

            This section of code processes mutually exclusive exon (MXE) events for a given gene. The goal is to identify cases where two 
            exons are mutually exclusive, meaning that one of the exons is included in some transcripts, while the other exon is included 
            in other transcripts of the same gene. These exons are never found together in the same transcript.

            The code iterates through all possible pairs of upstream and downstream exons, checking if the current exon is present in the 
            current transcript. If the current exon is not present, and there is an edge connecting the first and second flanking exons 
            with a candidate exon in between, a mutually exclusive event is detected.
            """

            # Iterate through the upstream exons
            for i in range(0, len(upstream)):
                # Extract the first flanking exon
                first_flanking = upstream[i]

                # Iterate through the downstream exons
                for j in range(0, len(downstream)):
                    # Extract the second flanking exon
                    second_flanking = downstream[j]

                    # Iterate through the transcripts in the gene
                    for transcript_id in genes[gene_id]:
                        # If we have the current exon, but not the first or second flanking exons, skip it, since this disqualifies it from being a mutually exclusive event
                        if [current_exon[0], current_exon[1]] in genes[gene_id][transcript_id] or first_flanking not in genes[gene_id][transcript_id] or second_flanking not in genes[gene_id][transcript_id]:
                            continue

                        # Extract the indices of the first and second flanking exons
                        first_flanking_index = genes[gene_id][transcript_id].index(first_flanking)
                        second_flanking_index = genes[gene_id][transcript_id].index(second_flanking)

                        # Candidate exon for a mutually exclusive event
                        mxe = genes[gene_id][transcript_id][first_flanking_index + 1]

                        # Check if the first and second flanking exons are adjacent and the candidate exon is not overlapping with the current exon, then it is a mutually exclusive event
                        if (first_flanking_index + 1 == second_flanking_index - 1) and (mxe[0] > current_exon[1]):
                            # Create a unique key for the mutually exclusive event
                            key = f"{sup_info[1]}:{current_exon[0] - 1}:{current_exon[1]}:{mxe[0] - 1}:{mxe[1]}:{first_flanking[1]}:{second_flanking[0] - 1}"

                            # Check if the key is already in the mutually exclusive events dictionary
                            if key in mx_events:
                                duplicate_mxe += 1

                            # Otherwise, add the key to the mutually exclusive events dictionary and increment the mutually exclusive event count
                            else:
                                num_mxe_events += 1
                                mx_events[key] = 1
                                mxe_output.write(f"{sup_info[1]}\t{gene_id}\t{mxe[0] - 1}\t{mxe[1]}\t{first_flanking[0] - 1}\t{first_flanking[1]}\t{sup_info[2]}\n")

            """
            2. ALTERNATIVE 3' (ALT-3) AND ALTERNATIVE 5' (ALT-5) SPLICE SITE EVENTS

            This section of code processes both alternative 3' (ALT-3) and alternative 5' (ALT-5) splice site events for a given gene. 
            The goal is to identify cases where there are alternative splice sites either at the 3' end or the 5' end of an exon. 
            This is done by examining the exons upstream and downstream flanking the exon of interest.

            The code iterates through all possible pairs of upstream and downstream exons, checking if they have the same start or end 
            positions but different end or start positions, indicating an alternative splice site event. If the strand is positive, it 
            differentiates between ALT-3 and ALT-5 events based on the start and end positions.

            For ALT-3 events, it checks pairs of downstream exons that share the same end position.
            For ALT-5 events, it checks pairs of upstream exons that share the same start position.
            """

            # Iterate through the upstream exons to find alternative splice sites
            for i in range(0, len(upstream) - 1):
                # Extract the first upstream exon
                first_exon = upstream[i]

                # Compare with other upstream exons
                for j in range(i + 1, len(upstream)):
                    # Extract the second upstream exon
                    second_exon = upstream[j]

                    # Check if the start positions are the same, indicating an alternative splice site event
                    if first_exon[0] == second_exon[0]:
                        # Create a unique key for the alternative splice site event
                        key = f"{sup_info[1]}:{min(first_exon[1], second_exon[1])}:{max(first_exon[1], second_exon[1])}:{current_exon[0] - 1}"

                        # Check if the strand is positive (indicating an ALT-5 event)
                        if sup_info[2] == "+":
                            # If the key is already in the ALT-5 events dictionary, increment the duplicate count and continue
                            if key in ss5_events:
                                duplicate_ss5 += 1
                            else:
                                # Otherwise, add the key to the ALT-5 events dictionary and increment the ALT-5 event count
                                ss5_events[key] = 1
                                num5 += 1
                                a5ss_output.write(f"{sup_info[1]}\t{gene_id}\t{first_exon[0] - 1}\t{min(first_exon[1], second_exon[1])}\t{current_exon[0] - 1}\t{current_exon[1]}\t{sup_info[2]}\n")
                        else:
                            # If the strand is negative (indicating an ALT-3 event)
                            # If the key is already in the ALT-3 events dictionary, increment the duplicate count and continue
                            if key in ss3_events:
                                duplicate_ss3 += 1
                            else:
                                # Otherwise, add the key to the ALT-3 events dictionary and increment the ALT-3 event count
                                ss3_events[key] = 1
                                num3 += 1
                                a3ss_output.write(f"{sup_info[1]}\t{gene_id}\t{first_exon[0] - 1}\t{min(first_exon[1], second_exon[1])}\t{current_exon[0] - 1}\t{current_exon[1]}\t{sup_info[2]}\n")

            # Iterate through the downstream exons to find alternative splice sites
            for i in range(0, len(downstream) - 1):
                # Extract the first downstream exon
                first_exon = downstream[i]

                # Compare with other downstream exons
                for j in range(i + 1, len(downstream)):
                    # Extract the second downstream exon
                    second_exon = downstream[j]

                    # Check if the end positions are the same, indicating an alternative splice site event
                    if first_exon[1] == second_exon[1]:
                        # Create a unique key for the alternative splice site event
                        key = f"{sup_info[1]}:{current_exon[1]}:{min(first_exon[0], second_exon[0]) - 1}:{max(first_exon[0], second_exon[0]) - 1}"

                        # Check if the strand is positive (indicating an ALT-3 event)
                        if sup_info[2] == "+":
                            # If the key is already in the ALT-3 events dictionary, increment the duplicate count and continue
                            if key in ss3_events:
                                duplicate_ss3 += 1
                            else:
                                # Otherwise, add the key to the ALT-3 events dictionary and increment the ALT-3 event count
                                ss3_events[key] = 1
                                num3 += 1
                                a3ss_output.write(f"{sup_info[1]}\t{gene_id}\t{max(first_exon[0], second_exon[0]) - 1}\t{first_exon[1]}\t{current_exon[0] - 1}\t{current_exon[1]}\t{sup_info[2]}\n")
                        else:
                            # If the strand is negative (indicating an ALT-5 event)
                            # If the key is already in the ALT-5 events dictionary, increment the duplicate count and continue
                            if key in ss5_events:
                                duplicate_ss5 += 1
                            else:
                                # Otherwise, add the key to the ALT-5 events dictionary and increment the ALT-5 event count
                                ss5_events[key] = 1
                                num5 += 1
                                a5ss_output.write(f"{sup_info[1]}\t{gene_id}\t{max(first_exon[0], second_exon[0]) - 1}\t{first_exon[1]}\t{current_exon[0] - 1}\t{current_exon[1]}\t{sup_info[2]}\n")

            """
            3. ALTERNATIVE FIRST EXON (AFE) AND ALTERNATIVE LAST EXON (ALE) EVENTS

            This section of code processes both alternative first exon (AFE) and alternative last exon (ALE) events for a given gene.
            The goal is to identify cases where there are alternative first or last exons in transcripts of the gene. 

            For AFE events, it checks if the current exon is the second exon in the transcript and compares the first exon with the 
            first exon of other transcripts to identify non-overlapping exons.
            For ALE events, it checks if the current exon is the second to last exon in the transcript and compares the last exon with 
            the last exon of other transcripts to identify non-overlapping exons.

            If the strand is positive, it differentiates between AFE and ALE events based on the position of the exons.
            """

            # This function determines if an exon is fully contained in an internal exon
            def fully_contained_in_internal_exon(myExon, myGeneID):
                # Iterate through the transcripts
                for mytranscript_id in genes[myGeneID]:
                    # For each internal exon
                    for intExon in genes[myGeneID][mytranscript_id][1:-1]:
                        # Check if the exon is fully contained
                        if intExon[0] <= myExon[0] and intExon[1] >= myExon[1]:  ## fully contained
                            return True
                return False

            # Iterate through the transcripts in the gene
            for transcript_id in genes[gene_id]:
                # Skip if there are less than two exons
                if len(genes[gene_id][transcript_id]) < 2:
                    continue

                # Check if the current exon is the second exon in the transcript (for AFE events)
                if [current_exon[0], current_exon[1]] == genes[gene_id][transcript_id][1]:
                    # Extract the first exon
                    first_exon = genes[gene_id][transcript_id][0]

                    # Skip if the first exon is fully contained in an internal exon
                    if fully_contained_in_internal_exon(first_exon, gene_id):
                        continue

                    # Iterate through the candidate transcripts
                    for candidate_transcript_id in genes[gene_id]:
                        # Skip if there are less than two exons in the candidate transcript
                        if len(genes[gene_id][candidate_transcript_id]) < 2:
                            continue

                        # Check if the current exon is the second exon in the candidate transcript
                        if [current_exon[0], current_exon[1]] == genes[gene_id][candidate_transcript_id][1]:
                            # Extract the candidate first exon
                            candidate_first_exon = genes[gene_id][candidate_transcript_id][0]

                            # Check if the first exon and candidate first exon are non-overlapping
                            if candidate_first_exon[0] > first_exon[1] or first_exon[0] > candidate_first_exon[1]:
                                # Check if the candidate first exon is fully contained in an internal exon
                                if fully_contained_in_internal_exon(candidate_first_exon, gene_id):
                                    continue

                                # Check if the strand is positive (indicating an AFE event)
                                if sup_info[2] == "+":
                                    # Create a unique key for the AFE event
                                    key = f"{sup_info[1]}:{current_exon[0] - 1}:{current_exon[1]}:{min(first_exon[0], candidate_first_exon[0]) - 1}:{min(first_exon[1], candidate_first_exon[1])}:{max(first_exon[0], candidate_first_exon[0]) - 1}:{max(first_exon[1], candidate_first_exon[1])}"

                                    # Check if the key is already in the AFE events dictionary
                                    if key in afe_events:
                                        pass  # Do nothing if it's a duplicate
                                    else:
                                        # Otherwise, add the key to the AFE events dictionary and increment the AFE event count
                                        afe_events[key] = 1
                                        num_afe += 1
                                        afe_output.write(f"{sup_info[1]}\t{gene_id}\t{max(first_exon[0], candidate_first_exon[0]) - 1}\t{max(first_exon[1], candidate_first_exon[1])}\t{current_exon[0] - 1}\t{current_exon[1]}\t{sup_info[2]}\n")
                                else:
                                    # If the strand is negative (indicating an ALE event)
                                    # Create a unique key for the ALE event
                                    key = f"{sup_info[1]}:{current_exon[0] - 1}:{current_exon[1]}:{min(first_exon[0], candidate_first_exon[0]) - 1}:{min(first_exon[1], candidate_first_exon[1])}:{max(first_exon[0], candidate_first_exon[0]) - 1}:{max(first_exon[1], candidate_first_exon[1])}"

                                    # Check if the key is already in the ALE events dictionary
                                    if key in ale_events:
                                        pass  # Do nothing if it's a duplicate
                                    else:
                                        # Otherwise, add the key to the ALE events dictionary and increment the ALE event count
                                        ale_events[key] = 1
                                        num_ale += 1
                                        ale_output.write(f"{sup_info[1]}\t{gene_id}\t{max(first_exon[0], candidate_first_exon[0]) - 1}\t{max(first_exon[1], candidate_first_exon[1])}\t{current_exon[0] - 1}\t{current_exon[1]}\t{sup_info[2]}\n")

                # Check if the current exon is the second to last exon in the transcript (for ALE events)
                if [current_exon[0], current_exon[1]] == genes[gene_id][transcript_id][-2]:
                    # Extract the last exon
                    last_exon = genes[gene_id][transcript_id][-1]

                    # Check if the last exon is fully contained in an internal exon
                    if fully_contained_in_internal_exon(last_exon, gene_id):
                        continue

                    # Iterate through the candidate transcripts
                    for candidate_transcript_id in genes[gene_id]:
                        # Skip if there are less than two exons in the candidate transcript
                        if len(genes[gene_id][candidate_transcript_id]) < 2:
                            continue

                        # Check if the current exon is the second to last exon in the candidate transcript
                        if [current_exon[0], current_exon[1]] == genes[gene_id][candidate_transcript_id][-2]:
                            # Extract the candidate last exon
                            candidate_last_exon = genes[gene_id][candidate_transcript_id][-1]

                            # Check if the last exon and candidate last exon are non-overlapping
                            if candidate_last_exon[0] > last_exon[1] or last_exon[0] > candidate_last_exon[1]:
                                # Check if the candidate last exon is fully contained in an internal exon
                                if fully_contained_in_internal_exon(candidate_last_exon, gene_id):
                                    continue

                                # Check if the strand is negative (indicating an AFE event)
                                if sup_info[2] == "-":
                                    # Create a unique key for the AFE event
                                    key = f"{sup_info[1]}:{current_exon[0] - 1}:{current_exon[1]}:{min(last_exon[0], candidate_last_exon[0]) - 1}:{min(last_exon[1], candidate_last_exon[1])}:{max(last_exon[0], candidate_last_exon[0]) - 1}:{max(last_exon[1], candidate_last_exon[1])}"

                                    # Check if the key is already in the AFE events dictionary
                                    if key in afe_events:
                                        pass  # Do nothing if it's a duplicate
                                    else:
                                        # Otherwise, add the key to the AFE events dictionary and increment the AFE event count
                                        afe_events[key] = 1
                                        num_afe += 1
                                        afe_output.write(f"{sup_info[1]}\t{gene_id}\t{min(last_exon[0], candidate_last_exon[0]) - 1}\t{min(last_exon[1], candidate_last_exon[1])}\t{current_exon[0] - 1}\t{current_exon[1]}\t{sup_info[2]}\n")
                                else:
                                    # If the strand is positive (indicating an ALE event)
                                    # Create a unique key for the ALE event
                                    key = f"{sup_info[1]}:{current_exon[0] - 1}:{current_exon[1]}:{min(last_exon[0], candidate_last_exon[0]) - 1}:{min(last_exon[1], candidate_last_exon[1])}:{max(last_exon[0], candidate_last_exon[0]) - 1}:{max(last_exon[1], candidate_last_exon[1])}"

                                    # Check if the key is already in the ALE events dictionary
                                    if key in ale_events:
                                        pass  # Do nothing if it's a duplicate
                                    else:
                                        # Otherwise, add the key to the ALE events dictionary and increment the ALE event count
                                        ale_events[key] = 1
                                        num_ale += 1
                                        ale_output.write(f"{sup_info[1]}\t{gene_id}\t{min(last_exon[0], candidate_last_exon[0]) - 1}\t{min(last_exon[1], candidate_last_exon[1])}\t{current_exon[0] - 1}\t{current_exon[1]}\t{sup_info[2]}\n")

            """
            4. RETAINED INTRON (RI) EVENTS

            This section of code processes retained intron (RI) events for a given gene. The goal is to identify cases where an intron is retained
            in some transcripts of the gene. 

            For RI events, it checks if there is an exon that spans from the upstream flanking exon to the end of the current exon, or 
            from the start of the current exon to the downstream flanking exon. If such an exon exists, it indicates a retained intron event.
            """

            # Iterate through the upstream exons to find retained intron events
            for i in range(0, len(upstream)):
                # Extract the first flanking exon
                first_flanking_exon = upstream[i]

                # Iterate through the transcripts in the gene
                for transcript_id in genes[gene_id]:
                    # Check if there is an exon that starts from the first flanking exon and ends at the current exon
                    if [first_flanking_exon[0], current_exon[1]] in genes[gene_id][transcript_id]:
                        # Create a unique key for the retained intron event
                        key = f"{sup_info[1]}:{first_flanking_exon[1]}:{current_exon[0] - 1}"

                        # Check if the key is already in the retained intron events dictionary
                        if key in ri_events:
                            duplicate_ri += 1
                            continue  # Move to the next transcript

                        # Otherwise, add the key to the retained intron events dictionary and increment the retained intron event count
                        else:
                            ri_events[key] = 1
                            num_ri += 1
                            ri_output.write(f"{sup_info[1]}\t{gene_id}\t{first_flanking_exon[0] - 1}\t{first_flanking_exon[1]}\t{sup_info[2]}\n")

            # Iterate through the downstream exons to find retained intron events
            for i in range(0, len(downstream)):
                # Extract the first flanking exon
                first_flanking_exon = downstream[i]

                # Iterate through the transcripts in the gene
                for transcript_id in genes[gene_id]:
                    # Check if there is an exon that starts from the current exon and ends at the first flanking exon
                    if [current_exon[0], first_flanking_exon[1]] in genes[gene_id][transcript_id]:
                        # Create a unique key for the retained intron event
                        key = sup_info[1] + ":" + str(current_exon[1]) + ":" + str(first_flanking_exon[0] - 1)

                        # Check if the key is already in the retained intron events dictionary
                        if key in ri_events:
                            duplicate_ri += 1
                            continue  # Move to the next transcript

                        # Otherwise, add the key to the retained intron events dictionary and increment the retained intron event count
                        else:
                            ri_events[key] = 1
                            num_ri += 1
                            ri_output.write(f"{sup_info[1]}\t{gene_id}\t{first_flanking_exon[0] - 1}\t{first_flanking_exon[1]}\t{sup_info[2]}\n")

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
    print_and_log("----------------------------------------------------------------", logger)
    print()
    logger.handlers[0].stream.write("\n\n")
