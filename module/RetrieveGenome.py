import os
import subprocess
import re
from module import JsonHelper
from Bio import SeqIO

STRAIN_FILE = "data/strains.json"

def getBarcode(filename):
    regex = re.compile('ESC_.+_AS')
    results = regex.findall(filename)
    num_of_results = len(results)
    if num_of_results == 1:
        return results[0]
    raise Warning("{} barcode is found".format(num_of_results))
    return ""

def getSerotype(barcode, strains):
    for strain in strains:
        if strain["assembly_barcode"] == barcode:
            return strain["serotype"]
    raise Warning("no serotype found by HeaderModifier")
    return ""

def getNewHeader(filename, strains):
    barcode = getBarcode(filename)
    serotype = getSerotype(barcode, strains)
    return barcode + '|' + serotype

def extend_fasta_header(filename, output_dir, strains):
    dst = output_dir + "/" + filename

    new_header = getNewHeader(filename, strains)
    new_records = []
    for record in SeqIO.parse(dst, "fasta"):
        record.id = new_header + "|" + record.id
        record.description = '' # remove duplicate description
        new_records.append(record)
    SeqIO.write(new_records, dst, "fasta")

def retrieve_genome(db_dir, input_dir, output_dir):
    strains = JsonHelper.read_from_json(STRAIN_FILE)
    genome_filenames_filtered = JsonHelper.read_from_json(input_dir)
    filtered_genome_filenames_len = len(genome_filenames_filtered)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    for index, filename in enumerate(genome_filenames_filtered):
        target_file_path = db_dir + filename
        copied_file_path = output_dir + filename
        cmd = str("cp " + target_file_path + " " + copied_file_path)
        subprocess.call(cmd, shell=True)
        extend_fasta_header(filename, output_dir, strains)
        print("{}/{} files copied and modified".format(index+1, filtered_genome_filenames_len))