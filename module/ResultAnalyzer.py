import logging
import re
from collections import defaultdict

import Bio
from Bio import SeqIO
from Bio.Blast import NCBIXML

import definitions
from module import JsonHelper, SerotypeHelper
import os
import logging

log = logging.getLogger(__name__)

BLACKLIST_FILE = os.path.join(definitions.OUTPUT_DIR, 'blacklist_genomes.json')
DICT_FILE = "output/genome_dict.json"
GENE_PAIRS = [('wzx', 'wzy'), ('wzt', 'wzm')]

def filterGenePair_helper(allele_dict, gene):
    for serotype, alleles in allele_dict.items():
        for allele in alleles:
            if allele['gene'] == gene:
                allele_dict[serotype].remove(allele)
    return allele_dict

def filterGenePair(allele_dict):
    '''
    Give a dict with structure:
        {serotype1:
            [{seq, gene},
            {seq, gene},
            {seq, gene},
            ...]
        {serotype1:
            [{seq, gene},
            {seq, gene},
            {seq, gene},
            ...]
        ...
        }
    if singles of gene pairs exist, filter them out.
    Then return the filtered allele_dict
    '''
    # get list of genes
    alleles = list(allele_dict.values())
    gene_list = [allele['gene']
                    for alleles in allele_dict.values()
                        for allele in alleles]
    for a, b in GENE_PAIRS:
        if (a in gene_list) and (b not in gene_list):
            allele_dict = filterGenePair_helper(allele_dict, a)
        if (b in gene_list) and (a not in gene_list):
            allele_dict = filterGenePair_helper(allele_dict, b)
    return allele_dict

def getGeneName(allele_desc):
    for gene_name in definitions.GENE_LIST:
        if gene_name in allele_desc:
            return gene_name
    log.info("No gene name found for %s", allele_desc)
    return ""

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

def initialize_serotype_dict():
    serotype_dict_file = 'output/serotype_dict.json'
    serotype_dict = defaultdict(list)
    allele_seqs = list(SeqIO.parse(definitions.SEROTYPED_ALLELE, 'fasta'))
    for allele_seq in allele_seqs:
        desc = allele_seq.description
        serotype = getSerotype(desc)
        if serotype == "":
            # These alleles have novel serotype. Ignored for now.
            continue
        new_entry = {
            "des": desc,
            "seq": str(allele_seq.seq),
            "gene": getGeneName(desc),
            "num": len(serotype_dict[serotype])+1
        }
        serotype_dict[serotype].append(new_entry)
    allele_count = sum(len(x) for x in serotype_dict.values())
    log.info("Serotype dictionary contains %d entries", allele_count)
    JsonHelper.write_to_json(serotype_dict, serotype_dict_file)
    return serotype_dict_file




def create_genome_result(blast_result_file):
    '''
    Return genome_output filepath
    '''
    genome_dict = defaultdict(list)
    blacklist_genomes = []
    blast_result = NCBIXML.parse(open(blast_result_file))
    if not blast_result:
        return
    for iteration in blast_result:
        if not iteration.alignments:
            continue
        allele_desc = iteration.query
        allele_serotype = getSerotype(allele_desc)
        for alignment in iteration.alignments:
            if not alignment.hsps:
                continue
            genome_desc = alignment.hit_def
            if genome_desc in blacklist_genomes:
                continue
            genome_serotypes = getSerotypes(genome_desc)
            is_mismatch = isMismatch(allele_serotype, genome_serotypes)
            is_same_class = isSameClass(allele_serotype, genome_serotypes)
            for hsp in alignment.hsps:
                match_len = hsp.identities
                query_len = iteration.query_length
                percent_identity = match_len / query_len
                if percent_identity < 0.97 or not is_same_class:
                    continue
                if match_len == query_len:
                    if is_mismatch:
                        log.debug("{0} and {1} are mismatch. {1} added to blacklist".format(allele_desc, genome_desc))
                        blacklist_genomes.append(genome_desc)
                    continue
                new_entry = {
                    "allele_desc": allele_desc,
                    "genome_desc": genome_desc,
                    "identity": "{}/{}".format(hsp.identities, iteration.query_length),
                    "Align Seq": {
                        "query": hsp.query,
                        "match": hsp.match,
                        "sbjct": hsp.sbjct
                    }
                }
                genome_dict[genome_desc].append(new_entry)
    # sort genome_dict by identity
    genome_dict_sorted = {}
    for genome_desc, alignments in genome_dict.items():
        #  filter out alleles associated with blacklisted genome
        if genome_desc in blacklist_genomes:
            continue
        genome_dict_sorted[genome_desc] = sorted(alignments, key=lambda alignment: alignment["identity"], reverse=True)
    genome_output_file = os.path.join(definitions.OUTPUT_DIR, 'genome_output.json')
    
    JsonHelper.write_to_json(genome_dict_sorted, genome_output_file)
    JsonHelper.write_to_json(blacklist_genomes, BLACKLIST_FILE)
    return genome_output_file


def expand_serotype_dict(serotype_dict_file, genome_dict_file):
    '''
    Expand serotype_dict.json by adding new alleles from blast result
    '''
    genome_dict = JsonHelper.read_from_json(genome_dict_file)
    serotype_dict = JsonHelper.read_from_json(serotype_dict_file)
    seq_hash = {}
    log.info("Start expanding serotype dictionary")
    # fill serotype_dict
    for genome_desc, alignments in genome_dict.items():
        new_dict = defaultdict(list)
        for alignment in alignments:
            new_item = {
                'seq':'',
                'gene':'',
                'des': "part of "+genome_desc.split("|")[0]}
            sbject = alignment['Align Seq']['sbjct']
            allele_desc = alignment["allele_desc"]
            serotype = getSerotype(allele_desc)
            # list of existing seqs in new_dict
            seq_list = [allele['seq']
                            for alleles in new_dict.values()
                                for allele in alleles]
            new_item['gene'] = getGeneName(allele_desc)
            # skip repeated seq
            if sbject in seq_list:
                continue
            new_item['seq'] = sbject
            new_dict[serotype].append(new_item)
        # make sure gene pair exist, otherwise, remove the singles
        filterGenePair(new_dict)
        # merge filtered dictionary to final dictionary
        for serotype, alleles in new_dict.items():
            for allele in alleles:
                # add only non-repeated sequence
                seq_list = [s_allele['seq']
                                for s_alleles in serotype_dict.values()
                                    for s_allele in s_alleles]
                new_seq = allele['seq']
                if new_seq not in seq_hash:
                    num = len(serotype_dict[serotype])+1
                    seq_hash[new_seq] = True
                    new_entry = {
                        'seq':allele['seq'],
                        'num':num,
                        'des':allele['des']}
                    serotype_dict[serotype].append(new_entry)

    allele_count = sum(len(x) for x in serotype_dict.values())
    log.info("Serotype dictionary contains %d entries", allele_count)
    JsonHelper.write_to_json(serotype_dict, "output/serotype_dict.json")
