import os
import subprocess
import re
from module import JsonHelper
import definitions
from Bio import SeqIO
import logging

log = logging.getLogger(__name__)

def getBarcode(filename):
    filename_no_ext = os.path.splitext(filename)[0]
    regex = re.compile('ESC_.+_AS')
    results = regex.findall(filename_no_ext)
    num_of_results = len(results)
    if num_of_results == 1:
        return results[0]
    log.debug("no serotype found by getBarcode()")
    return ""

def getSerotype(barcode, strains):
    for strain in strains:
        if strain["assembly_barcode"] == barcode:
            return strain["serotype"]
    log.debug("no serotype found by getSerotype()")
    return ""

def getNewHeader(filename, strains):
    barcode = getBarcode(filename)
    serotype = getSerotype(barcode, strains)
    return '|'.join([barcode, serotype])

def extend_fasta_header(filepath, strains):
    filename = os.path.basename(filepath)
    new_header = getNewHeader(filename, strains)
    new_records = []
    for record in SeqIO.parse(filepath, "fasta"):
        record.id = new_header + "|" + record.id
        record.description = '' # remove duplicate description
        new_records.append(record)
    SeqIO.write(new_records, filepath, "fasta")

def retrieve_genomes(strain_file, selected_genomes_json):
    strains = JsonHelper.read_from_json(strain_file)
    selected_genomes = JsonHelper.read_from_json(selected_genomes_json)
    genome_num = len(selected_genomes)
    if not os.path.exists(definitions.GENOME_DIR):
        os.makedirs(definitions.GENOME_DIR)
    for index, filename in enumerate(selected_genomes):
        target_filepath = os.path.join(definitions.DB_DIR, filename)
        copied_filepath = os.path.join(definitions.GENOME_DIR, filename)
        cmd = ["cp", target_filepath, copied_filepath]
        subprocess.run(cmd, check=True)
        extend_fasta_header(copied_filepath, strains)
        print("{}/{} files copied and modified".format(index+1, genome_num))