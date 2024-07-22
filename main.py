from methods.preprocess import preprocess_gene_ids
from methods.as_events import find_as_events
import configparser
import logging
import os
import datetime
from methods.print_and_log import print_and_log

# Read in the configuration values
config = configparser.ConfigParser()
config.read("./config/config.ini")
gtf_file = config["INPUT_FILES"]["gtf_file"]
ref_file = config["INPUT_FILES"]["ref_file"]
as_output_dir = config["OUTPUT"]["as_output_dir"]
gtf_output_file = config["OUTPUT"]["gtf_output_file"]
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

# Preprocess the Gene IDs, creating a correct GTF file
preprocess_gene_ids(gtf_file, ref_file, gtf_output_file, logger)

# Find alternative splicing events with the new GTF file
find_as_events(gtf_output_file, gtf_output_file.split("/")[-1].split(".")[0], as_output_dir, logger)
