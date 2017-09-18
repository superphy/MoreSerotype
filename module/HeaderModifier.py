import os
import re
import json
import subprocess
from Bio import SeqIO
from module import JsonHelper

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

def extend_fasta_header(filename, strains, genome_dir, output_dir):
    src = genome_dir + "/" + filename
    dst = output_dir + "/" + filename

    new_header = getNewHeader(filename, strains)
    cmd = str("cp " + src + " " + dst)
    subprocess.call(cmd, shell=True)
    new_records = []
    for record in SeqIO.parse(dst, "fasta"):
        record.id = new_header + "|" + record.id
        record.description = '' # remove duplicate description
        new_records.append(record)
    SeqIO.write(new_records, dst, "fasta")
    print("new file created at " + dst)

def modify_header(genome_dir, output_dir):
    strains = JsonHelper.read_from_json(STRAIN_FILE)
    genome_filenames = os.listdir(genome_dir)
    genome_filenames_num = len(genome_filenames)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    for index, filename in enumerate(genome_filenames):
        print("At iteration number {}/{}". format(index+1, genome_filenames_num))
        extend_fasta_header(filename, strains, genome_dir, output_dir)