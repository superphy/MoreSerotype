from module import JsonHelper
import os

STRAIN_FILE = "data/strains.json"

def filter_genome(genome_dir, filtered_file):
    print("Start filtering genome")
    strains = JsonHelper.read_from_json(STRAINS_FILE)
    genome_filenames_filtered = []
    genome_filenames = os.listdir(genome_dir)
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
    JsonHelper.write_to_json(genome_filenames_filtered, filtered_file)