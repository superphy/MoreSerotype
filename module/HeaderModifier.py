import os
import re
import json
import subprocess
from Bio import SeqIO
from module import JsonHelper
INPUT_DIR = "data/filtered_genome"
OUTPUT_DIR = "data/modified_genome"
INPUT_FILE = "data/strains.json"

def getBarcode(filename):
    regex = re.compile('ESC_.+_AS')
    results = regex.findall(filename)
    num_of_results = len(results)
    if num_of_results == 1:
        return results[0]
    else:
        raise Exception("{} barcode is found".format(num_of_results))
    return ""

def getSerotype(barcode, strains):
    for strain in strains:
        if strain["assembly_barcode"] == barcode:
            return strain["serotype"]
    print("no serotype found by HeaderModifier")    
    return ""

def extend_fasta_header(filename, strains):
    src = INPUT_DIR + "/" + filename
    dst = OUTPUT_DIR + "/" + filename

    barcode = getBarcode(filename)
    serotype = getSerotype(barcode, strains)

    new_header = barcode + '|' + serotype
    cmd = str("cp " + src + " " + dst)
    subprocess.call(cmd, shell=True)
    new_records = []
    for record in SeqIO.parse(dst, "fasta"):
        record.id = new_header + "|" + record.id
        record.description = '' # remove duplicate description
        new_records.append(record)
    SeqIO.write(new_records, dst, "fasta")
    print("new file created at " + dst)

def modify_header():
    strains = JsonHelper.read_from_json(INPUT_FILE)
    genome_filenames = os.listdir(INPUT_DIR)
    genome_filenames_num = len(genome_filenames)
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    for index, filename in enumerate(genome_filenames):
        print("At iteration number {}/{}". format(index+1, genome_filenames_num))
        extend_fasta_header(filename, strains)

