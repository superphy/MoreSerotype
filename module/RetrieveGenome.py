from module import JsonHelper
import os
import subprocess

def retrieve_genome(db_dir, input_dir, output_dir):
    genome_filenames_filtered = JsonHelper.read_from_json(input_dir)
    filtered_genome_filenames_len = len(genome_filenames_filtered)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    for index, filename in enumerate(genome_filenames_filtered):
        target_file_path = db_dir + filename
        copied_file_path = output_dir + filename
        cmd = str("cp " + target_file_path + " " + copied_file_path)
        subprocess.call(cmd, shell=True)
        print("{}/{} files copied".format(index+1, filtered_genome_filenames_len))