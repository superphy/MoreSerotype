from collections import defaultdict
from module import JsonHelper
from Bio import SeqIO
import re

ALLELE_FILE = "data/EcOH.fasta"
DICT_FILE = "output/serotype_dict.json"

def getSerotypes(str):
    serotypes = {'O':'', 'H':''}
    regex = re.compile("(O\d{1,3})(?!\d)")
    results = regex.findall(str)
    results_len = len(results)
    if results_len > 0:
        serotypes['O'] = results[0][0]
    regex = re.compile("(H\d{1,3})(?!\d)")
    results = regex.findall(str)
    results_len = len(results)
    if results_len > 0:
        serotypes['H'] = results[0][0]
    return serotypes


def getSerotype(str):
    serotypes = {'O':'', 'H':''}
    regex = re.compile("((O|H)\d{1,3})(?!\d)")
    results = regex.findall(str)
    results_len = len(results)
    if results_len > 0:
        return results[0][0]
    return ''

def isMismatch(allele_serotype, genome_serotype):
    if not allele_serotype:
        return False
    serotype_type = allele_serotype[0]
    if genome_serotype[serotype_type] == '':
        return False
    if allele_serotype != serotype_type + genome_serotype[serotype_type]:
        print(allele_serotype, genome_serotype, " are mismatch")
        return True
    return False

def isSameClass(allele_serotype, genome_serotype):
    if not allele_serotype:
        return False
    serotype_type = allele_serotype[0]
    if genome_serotype[serotype_type] != '':
        return True
    return False

def initialize_dict(allele_file=ALLELE_FILE, dict_file=DICT_FILE):
    serotype_dict = defaultdict(list)
    allele_seqs = list(SeqIO.parse(allele_file, 'fasta'))
    for allele_seq in allele_seqs:
        desc = allele_seq.description
        serotype = getSerotype(desc)
        if serotype == "":
            continue
        new_entry = {
            "desc": allele_seq.description,
            "seq": str(allele_seq.seq)
        }
        serotype_dict[serotype].append(new_entry)
    JsonHelper.write_to_json(serotype_dict, dict_file)