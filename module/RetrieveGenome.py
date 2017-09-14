from module import JsonHelper
import os
import subprocess

DATABASE_PATH = '~/moria/enterobase_db/'
INPUT_FILE = "data/genome_filenames_filtered.json"
OUTPUT_DIR = "data/filtered_genome/"

def retrieve_genome():
    genome_filenames_filtered = JsonHelper.read_from_json(INPUT_FILE)
    filtered_genome_filenames_len = len(genome_filenames_filtered)
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    for index, filename in enumerate(genome_filenames_filtered):
        target_file_path = DATABASE_PATH + filename
        copied_file_path = OUTPUT_DIR + filename
        cmd = str("cp " + target_file_path + " " + copied_file_path)
        subprocess.call(cmd, shell=True)
        print("{}/{} files copied".format(index+1, filtered_genome_filenames_len))