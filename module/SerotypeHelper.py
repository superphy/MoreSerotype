from collections import defaultdict
from module import JsonHelper
from Bio import SeqIO
import re


def getSerotypes(str):
    serotypes = {'O':'', 'H':''}
    for key, value in serotypes.items():
        regex = re.compile("(?<!(Non-))("+key+"\d{1,3})(?!\d)")
        results = regex.findall(str)
        results_len = len(results)
        if results_len > 0:
            serotypes[key] = results[0][1][1:]
    return serotypes


def getSerotype(str):
    serotypes = ""
    regex = re.compile("(?<!(Non-))((O|H)\d{1,3})(?!\d)")
    results = regex.findall(str)
    results_len = len(results)
    if results_len > 0:
        return results[0][1]
    return ''

def isMismatch(allele_serotype, genome_serotype):
    if not allele_serotype:
        return False
    serotype_type = allele_serotype[0]
    if genome_serotype[serotype_type] == '':
        return False
    if allele_serotype != serotype_type + genome_serotype[serotype_type]:
        # print(allele_serotype, genome_serotype, " are mismatch")
        return True
    return False

def isSameClass(allele_serotype, genome_serotype):
    if not allele_serotype:
        return False
    serotype_type = allele_serotype[0]
    if genome_serotype[serotype_type] != '':
        return True
    return False

def initialize_dict(allele_file, dict_file):
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