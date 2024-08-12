from methods.preprocess import preprocess_gene_ids
from methods.preprocess import preprocess_bam_files
from methods.as_quant import as_quant
from methods.find_as_events import find_as_events
import configparser
import logging
import os
import datetime
from methods.print_and_log import print_and_log
import time

# Read in the configuration values
config = configparser.ConfigParser()
config.read("./config/config.ini")

# Input values
gtf_file = config["INPUT_FILES"]["gtf_file"]
ref_file = config["INPUT_FILES"]["ref_file"]
group1_folder = config["INPUT_FILES"]["group1_folder"]
group2_folder = config["INPUT_FILES"]["group2_folder"]

# Output values
as_events_output_dir = config["OUTPUT"]["as_events_output_dir"]
as_quant_output_dir = config["OUTPUT"]["as_quant_output_dir"]
annotation_file = config["OUTPUT"]["annotation_file"]
gtf_output_file = config["OUTPUT"]["gtf_output_file"]

# Parameters
samtools = config["PARAMETERS"]["samtools_dir"]
cores = config["PARAMETERS"]["cores"]
species = config["PARAMETERS"]["species"]
method = config["PARAMETERS"]["method"]
novel = config["PARAMETERS"]["novel"]


# Logging Folder
log_dir = config["LOGGING"]["log_dir"]

# Configure the logging object
log_file_name = os.path.join(log_dir, f"as_quant.{datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.log")
logging.basicConfig(filename=log_file_name, level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger()

# Begin the pipeline
print_and_log("----------------------------------------------------------------", logger)
print_and_log("| Beginning AS-Quant Pipeline...                               |", logger)
print_and_log("----------------------------------------------------------------", logger)
logger.handlers[0].stream.write("\n\n")
print()

# Print the config file for logging
print_and_log("----------------------------------------------------------------", logger)
print_and_log("| Configuration File:                                          |", logger)
print_and_log("----------------------------------------------------------------", logger)
for section in config:
    if section == "DEFAULT":
        continue
    print_and_log(section, logger)
    for key in config[section]:
        print_and_log(f"{key}: {config[section][key]}", logger)
    print()
    logger.handlers[0].stream.write("\n")
print_and_log("----------------------------------------------------------------", logger)
logger.handlers[0].stream.write("\n\n")
print()


# Track time
start = time.perf_counter()

# Preprocess the Gene IDs, creating a correct GTF file
prestart = time.perf_counter()
preprocess_gene_ids(gtf_file, ref_file, gtf_output_file, annotation_file, logger)
print_and_log(f"Step 0 Completed in {time.perf_counter() - prestart:.2f} seconds", logger)
print_and_log("----------------------------------------------------------------", logger)
logger.handlers[0].stream.write("\n\n")
print()

# Find alternative splicing events with the new GTF file
find_as_events(gtf_output_file, species, as_events_output_dir, logger)

# Preprocess the BAM files using samtools
preprocess_bam_files(group1_folder, group2_folder, species, samtools, cores, logger)

# Run the AS_Quant process
as_quant(species, as_quant_output_dir, as_events_output_dir, group1_folder, group2_folder, annotation_file, method, novel, cores, logger)

# Begin the pipeline
print_and_log("----------------------------------------------------------------", logger)
print_and_log(f"| AS-Quant pipeline completed in {time.perf_counter() - start:.2f} seconds                 |", logger)
print_and_log("----------------------------------------------------------------", logger)
