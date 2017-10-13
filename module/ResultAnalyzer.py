import logging
import os
import re
from collections import defaultdict

import Bio
from Bio import SeqIO
from Bio.Blast import NCBIXML

import definitions
from module import JsonHelper

LOG = logging.getLogger(__name__)

BLACKDICT_FILE = os.path.join(definitions.OUTPUT_DIR, 'blacklist_dict.json')
DICT_FILE = "output/genome_dict.json"
GENE_PAIRS = [('wzx', 'wzy'), ('wzt', 'wzm')]

def filter_lone_pair_helper(allele_dict, gene):
    '''
    Remove lone pair from the dictionary
    '''
    for serotype, alleles in allele_dict.items():
        for allele in alleles:
            if allele['gene'] == gene:
                allele_dict[serotype].remove(allele)
    return allele_dict

def filter_lone_pair(serotype_dict):
    '''
    {
        '1':[
            {'des':'','seq':'','gene':'wzy'}
        ]
    }
    -->
    {
        '1':{}
    }
    Args:
        serotype_dict(dict): serotype dictionary
    
    Returns:
        serotype_dict(dict): filteredserotype dictionary
    '''
    gene_list = [
        allele['gene']
            for alleles in list(serotype_dict.values())
                for allele in alleles
    ]
    for a, b in GENE_PAIRS:
        if a in gene_list and b in gene_list:
            continue
        elif a not in gene_list and b not in gene_list:
            continue
        else:
            LOG.debug("Mismatch occured. Remove %s and %s from current iteration.", a, b)
            serotype_dict = filter_lone_pair_helper(serotype_dict, a)
            serotype_dict = filter_lone_pair_helper(serotype_dict, b)
    return serotype_dict

def getGeneName(str):
    '''
    Args:
        str(string): Allele or genome description

    Return:
        string: Gene name if found, '' if not found
    '''
    for gene_name in definitions.GENE_LIST:
        if gene_name in str:
            return gene_name
    LOG.info("No gene name found for %s", str)
    return ""

def getSerotypes(str):
    '''
    Args:
        str(string): Allele or genome description

    Return:
        dict: {'O': [string], 'H': [string]}
    '''
    serotypes = {'O': '','H': ''}
    for key in ['O','H']:
        regex = re.compile("(?<!(Non-))("+key+"\d{1,3})(?!\d)")
        results = regex.findall(str)
        results_len = len(results)
        if results_len > 0:
            serotypes[key] = results[0][1][1:]
    return serotypes

def getAntigen(serotypes):
    '''
    Args:
        serotypes(dict): {'O': [string], 'H': [string]}
    Returns:
        (string): 'O' or 'H' or ''
    '''
    for antigen in ['O','H']:
        if serotypes[antigen] != '':
            return antigen
    return ''

def create_genome_result(blast_output_file):
    '''
    Return genome_output filepath
    '''
    LOG.info("Parsing blast results in %s", blast_output_file)

    genome_dict = defaultdict(list)
    blacklist_dict = {}
    result_handle = open(blast_output_file, 'r')

    for line in result_handle:
        fields = line.strip().split()

        #  '6 qseqid qlen qseq sseqid length sseq pident'
        blast_record = {
            'qseqid': fields[0],
            'qlen': fields[1],
            'qseq': fields[2],
            'sseqid': fields[3],
            'length': fields[4],
            'sseq': fields[5],
            'pident': fields[6]
        }

        allele_name = blast_record['qseqid']
        allele_serotypes = getSerotypes(blast_record['qseqid'])
        allele_antigen = getAntigen(allele_serotypes)
        genome_desc = blast_record['sseqid']
        genome_name = genome_desc.split("|")[0]
        genome_serotypes = getSerotypes(blast_record['sseqid'])
        identity = float(blast_record['pident']) / 100.0
        length = blast_record['length']
        qlen = blast_record['qlen']
        
        if allele_antigen == "":
            # These alleles have novel serotype. Ignored for now.
            continue
        if genome_serotypes[allele_antigen] == '':
            # allele and genome have different antigen
            continue
        if allele_serotypes[allele_antigen] != genome_serotypes[allele_antigen]:
            # allele and genome serotype mismatch
            if identity >= 1 and length==qlen:
                # perfect mismatch
                LOG.debug("{0} and {1} are mismatch. {1} added to blacklist".format(allele_name, genome_name))
                new_entry = {
                    'genome name': genome_name,
                    'allele name': allele_name,
                    'given O': genome_serotypes['O'],
                    'given H': genome_serotypes['H'],
                    'predicted O': allele_serotypes['O'],
                    'predicted H': allele_serotypes['H']
                }
                if genome_name in blacklist_dict:
                    # append to existing entry
                    existing_entry = blacklist_dict[genome_name]
                    new_entry = merge_blacklist_entry(existing_entry, new_entry)
                blacklist_dict[genome_name]=new_entry
                continue
        else:
            # matching serotype
            if identity < 1.0:
                # non-perfect match
                new_entry = {
                    "allele_name": allele_name,
                    "genome_name": genome_name,
                    "identity": identity,
                    "Align Seq": {
                        "query": blast_record['qseq'],
                        "sbjct": blast_record['sseq']
                    }
                }
                genome_dict[genome_name].append(new_entry)
    # sort genome_dict by identity
    genome_dict_sorted = {}
    for genome_desc, alignments in genome_dict.items():
        #  filter out alleles associated with blacklisted genome
        genome_dict_sorted[genome_desc] = sorted(alignments, key=lambda alignment: alignment["identity"], reverse=True)
    blacklist_list = list(blacklist_dict.values())
    genome_output_file = os.path.join(definitions.OUTPUT_DIR, 'genome_output.json')
    blacklist_dict_file = os.path.join(definitions.OUTPUT_DIR, 'blacklist_dict.json')
    JsonHelper.write_to_json(blacklist_list, blacklist_dict_file)
    JsonHelper.write_to_json(genome_dict_sorted, genome_output_file)
    return genome_output_file


def expand_serotype_dict(genome_dict_file):
    '''
    Expand serotype_dict.json by adding new alleles from blast result
    '''
    serotype_dict = initialize_serotype_dict()
    genome_dict = JsonHelper.read_from_json(genome_dict_file)
    seq_hash = {}
    LOG.info("Start expanding serotype dictionary")
    for alignments in genome_dict.values():
        # create sub-dictionary for this genome
        new_dict = defaultdict(list)
        for alignment in alignments:
            allele_name = alignment["allele_name"]
            genome_name = alignment['genome_name']
            sbject = alignment['Align Seq']['sbjct']
            gene_name = getGeneName(allele_name)
            serotypes = getSerotypes(allele_name)
            antigen = getAntigen(serotypes)
            serotype = serotypes[antigen]
            serotype_tag = antigen+serotype
            # skip repeated seq
            if sbject in seq_hash:
                continue
            else:
                seq_hash[sbject] = True
            new_entry = {
                'gene': gene_name,
                'seq': sbject,
                'name': "part_of_"+genome_name
            }
            new_dict[serotype_tag].append(new_entry)
        # make sure gene pair exist, otherwise, remove the singles
        new_dict = filter_lone_pair(new_dict)
        # merge filtered dictionary to final dictionary
        serotype_dict = merge_serotype_dict(new_dict, serotype_dict)
    allele_count = sum(len(x) for x in serotype_dict.values())
    LOG.info("Serotype dictionary contains %d entries", allele_count)
    JsonHelper.write_to_json(serotype_dict, "output/serotype_dict.json")
    return allele_count

def initialize_serotype_dict():
    '''
    initialize_serotype_dict
    '''
    serotype_dict = defaultdict(list)
    allele_seqs = list(SeqIO.parse(definitions.SEROTYPED_ALLELE, 'fasta'))
    for allele_seq in allele_seqs:
        allele_name = allele_seq.description
        serotypes = getSerotypes(allele_name)
        antigen = getAntigen(serotypes)
        if antigen == "":
            # These alleles have novel serotype. Ignored for now.
            continue
        serotype = serotypes[antigen]
        serotype_tag = antigen+serotype
        new_entry = {
            "name": allele_name,
            "seq": str(allele_seq.seq),
            "gene": getGeneName(allele_name),
            "num": len(serotype_dict[serotype_tag])+1
        }
        serotype_dict[serotype_tag].append(new_entry)
    allele_count = sum(len(x) for x in serotype_dict.values())
    LOG.info("Serotype dictionary contains %d entries", allele_count)
    return serotype_dict

def merge_serotype_dict(new_dict, old_dict):
    '''
    Merge serotype dictionary
    '''
    # Get all existing allele names
    seq_list = [s_allele['seq']
        for s_alleles in old_dict.values()
            for s_allele in s_alleles]
    for serotype, alleles in new_dict.items():
        for allele in alleles:
            # add only non-repeated sequence
            num = len(old_dict[serotype])+1
            new_entry = {
                'seq':allele['seq'],
                'gene':allele['gene'],
                'num':num,
                'name':allele['name']
            }
            if new_entry['seq'] in seq_list:
                continue
            else:
                seq_list.append(new_entry['seq'])
                old_dict[serotype].append(new_entry)
    return old_dict

def merge_blacklist_entry(existing_entry, new_entry):
    existing_allele = existing_entry['allele name']
    new_allele = new_entry['allele name']
    for antigen in ['O', 'H']:
        given_serotype = existing_entry['given '+antigen]
        new_serotype = new_entry['predicted '+antigen]
        if new_serotype is '':
            continue
        old_serotype = existing_entry['predicted '+antigen]
        if new_serotype in old_serotype.split('|'):
            continue
        if old_serotype is '':
            existing_entry['predicted '+antigen] = new_serotype
            existing_entry['allele name'] = '|'.join([existing_allele, new_allele])
            continue
        # multiple prediction
        existing_entry['predicted '+antigen] = '|'.join([old_serotype, new_serotype])
        existing_entry['allele name'] = '|'.join([existing_allele, new_allele])
        LOG.warning("multiple predictions found for genome,%s", existing_entry['genome name'])
        LOG.warning("%s antigen, old: %s %s, new: %s, given: %s", antigen, existing_allele, old_serotype, new_serotype, given_serotype)
    return existing_entry