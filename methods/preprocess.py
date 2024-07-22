import sys
import time
from methods.print_and_log import print_and_log


def preprocess_gene_ids(gtf, ref, outf, logger):
    # Step 0: Correct gene_id in GTF file using refFlat file
    print_and_log("----------------------------------------------------------------", logger)
    print_and_log("| STEP 0: Correcting gene_id in GTF file using refFlat file... |", logger)
    print_and_log("----------------------------------------------------------------", logger)

    start = time.perf_counter()

    # Create a dictionary to store gene_id to gene_name mappings from refFlat file
    hash_dict = {}
    with open(ref, "r") as ref_file:
        for line in ref_file:
            info = line.split()
            # info[1] is gene_id and info[0] is gene_name in refFlat file
            hash_dict[info[1]] = info[0].replace('"', "").replace(";", "")
    # print(hash_dict)

    # Create a dictionary to record transcripts that don't have gene ID in refFlat
    list_dict = {}

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
                # print("Gene ID Extracted:", g_id)

                # Check for XM, XR, and YP prefixes in gene_id, and filter them out
                if g_id.startswith("XM") or g_id.startswith("XR") or g_id.startswith("YP"):
                    drop_count += 1
                    continue

                # If gene_id is in refFlat hash dictionary, use the corresponding gene_name
                if g_id in hash_dict:
                    # print("Gene ID Found in refFlat, correcting:", gene_id)
                    gene_id = hash_dict[g_id]
                else:
                    # If gene_id is not in refFlat, drop the transcript
                    drop_count += 1
                    continue

                # Update the attributes field with the new or original gene_id
                info[8] = f'gene_id "{gene_id}"; transcript_id "{g_id}";'
            else:
                print(f"Wrong line in GTF file:\n{line}", file=sys.stderr)
                sys.exit(1)

            # Write the updated line to the output file
            count += 1
            out_file.write("\t".join(info) + "\n")

    # Print the number of transcripts without a gene ID to stderr
    n = len(list_dict)
    print_and_log(f"{count} transcripts processed ", logger)
    print_and_log(f"{drop_count} transcripts dropped ", logger)
    print_and_log(f"Step 0 Completed in {time.perf_counter() - start:.2f} seconds", logger)
    print_and_log("----------------------------------------------------------------", logger)
    logger.handlers[0].stream.write('\n\n')

    # # Print each missing gene_id to stderr
    # for key in list_dict:
    #     print(f"~{key}~", file=sys.stderr)
