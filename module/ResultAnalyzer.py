from collections import defaultdict
from Bio.Blast import NCBIXML
from module import SerotypeHelper, JsonHelper
import logging

BLACKLIST_FILE = "output/blacklist_genomes.json"
DICT_FILE = "output/genome_dict.json"

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

def addUsefulAllele():
    