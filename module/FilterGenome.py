from module import JsonHelper
import os

INPUT_FILE = 'data/strains.json'
OUTPUT_FILE = 'data/genome_filenames_filtered.json'
GENOME_FILES_DIR = '/home/sam/moria/enterobase_db'




def filter_genome():
    print("Start filtering genome")
    strains = JsonHelper.read_from_json(INPUT_FILE)
    genome_filenames_filtered = []
    genome_filenames = os.listdir(GENOME_FILES_DIR)
    genome_filenames_len = len(genome_filenames)
    for index, filename in enumerate(genome_filenames):
        print("{}/{} genome files filtered".format(index, genome_filenames_len))
        for strain in strains:
            serotype = strain["serotype"]
            assembly_barcode = strain['assembly_barcode']
            if serotype == '' or assembly_barcode == '':
                print("something is wrong")
                continue
            if assembly_barcode in filename:
                genome_filenames_filtered.append(filename)
    JsonHelper.write_to_json(genome_filenames_filtered, OUTPUT_FILE)