from methods.preprocess import preprocess_gene_ids
from methods.as_events import find_as_events
import configparser

# Read in the configuration values
config = configparser.ConfigParser()
config.read("./config/config.ini")
gtf_file = config["INPUT_FILES"]["gtf_file"]
ref_file = config["INPUT_FILES"]["ref_file"]
as_output_dir = config["OUTPUT"]["as_output_dir"]
gtf_output_file = config["OUTPUT"]["gtf_output_file"]
log_dir = config["LOGGING"]["log_dir"]

# Preprocess the Gene IDs, creating a correct GTF file
preprocess_gene_ids(gtf_file, ref_file, gtf_output_file)

# Find alternative splicing events with the new GTF file
find_as_events(gtf_output_file, gtf_output_file.split("/")[-1].split(".")[0], as_output_dir, log_dir)
