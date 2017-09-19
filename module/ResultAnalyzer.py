from collections import defaultdict
from Bio.Blast import NCBIXML
from module import SerotypeHelper, JsonHelper
import logging
import Bio

BLACKLIST_FILE = "output/blacklist_genomes.json"
DICT_FILE = "output/genome_dict.json"
GENE_PAIRS = [('wzx','wzy'),('wzt','wzm')]

def create_genome_result(result_file):
    genome_dict = defaultdict(list)
    blacklist_genomes = []
    blast_result = NCBIXML.parse(open(result_file))
    if not blast_result:
        return
    for iteration in blast_result:
        if not iteration.alignments:
            continue
        allele_desc = iteration.query
        allele_serotype = SerotypeHelper.getSerotype(allele_desc)
        for alignment in iteration.alignments:
            if not alignment.hsps:
                continue
            genome_desc = alignment.hit_def
            if genome_desc in blacklist_genomes:
                continue
            genome_serotypes = SerotypeHelper.getSerotypes(genome_desc)
            is_mismatch = SerotypeHelper.isMismatch(allele_serotype, genome_serotypes)
            is_same_class = SerotypeHelper.isSameClass(allele_serotype, genome_serotypes)
            for hsp in alignment.hsps:
                match_len = hsp.identities
                query_len = iteration.query_length
                percent_identity = match_len / query_len
                if percent_identity < 0.97 or not is_same_class:
                    continue
                if match_len == query_len:
                    if is_mismatch:
                        print("{0} and {1} are mismatch. {1} added to blacklist".format(allele_desc, genome_desc))
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
        
    JsonHelper.write_to_json(genome_dict_sorted,DICT_FILE)
    JsonHelper.write_to_json(blacklist_genomes, BLACKLIST_FILE)

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
    for a, b in GENE_PAIRS:
        if a in allele_desc:
            return a
        if b in allele_desc:
            return b
    return ""

def addUsefulAllele():
    genome_dict = JsonHelper.read_from_json(DICT_FILE)
    serotype_dict = JsonHelper.read_from_json("output/serotype_dict.json")
    seq_hash = {}
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
            serotype = SerotypeHelper.getSerotype(allele_desc)
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
    JsonHelper.write_to_json(serotype_dict, "output/serotype_dict.json")