from methods.preprocess import preprocess_gene_ids
import configparser

# Read in the configuration values
config = configparser.ConfigParser()
config.read("./config/config.ini")
gtf_file = config["INPUT_FILES"]["gtf_file"]
ref_file = config["INPUT_FILES"]["ref_file"]
output_file = config["OUTPUT"]["output_file"]

preprocess_gene_ids(gtf_file, ref_file, output_file)
