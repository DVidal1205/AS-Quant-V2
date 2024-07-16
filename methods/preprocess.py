import sys


def preprocess_gene_ids(gtf, ref, outf):
    # Create a dictionary to store gene_id to gene_name mappings from refFlat file
    hash_dict = {}
    with open(ref, "r") as ref_file:
        for line in ref_file:
            info = line.split()
            # info[1] is gene_id and info[12] is gene_name in refFlat file
            hash_dict[info[0]] = info[11]

    # Create a dictionary to record transcripts that don't have gene ID in refFlat
    list_dict = {}

    # Open the GTF file for reading and the output file for writing
    with open(gtf, "r") as gtf_file, open(outf, "w") as out_file:
        for line in gtf_file:
            line = line.strip()  # Remove leading/trailing whitespace
            info = line.split("\t")  # Split line by tab to get GTF fields

            # Check if the line contains gene_id and transcript_id
            if 'gene_id "' in info[8] and 'transcript_id "' in info[8]:
                gene_id = ""
                # Extract gene_id and transcript_id from the attributes field (info[8])
                g_id = info[8].split('gene_id "')[1].split('";')[0]

                # If gene_id is in refFlat hash dictionary, use the corresponding gene_name
                if g_id in hash_dict:
                    gene_id = hash_dict[g_id]
                else:
                    # If gene_id is not in refFlat, keep the original gene_id
                    gene_id = g_id
                    list_dict[g_id] = 1  # Record the missing gene_id

                # Update the attributes field with the new or original gene_id
                info[8] = f'gene_id "{gene_id}"; transcript_id "{g_id}";'
            else:
                print(f"Wrong line in GTF file:\n{line}", file=sys.stderr)
                sys.exit(1)

            # Write the updated line to the output file
            out_file.write("\t".join(info) + "\n")

    # Print the number of transcripts without a gene ID to stderr
    n = len(list_dict)
    print(f"There are {n} transcripts that don't have gene id.", file=sys.stderr)

    # Print each missing gene_id to stderr
    for key in list_dict:
        print(f"~{key}~", file=sys.stderr)
