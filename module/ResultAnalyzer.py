from collections import defaultdict
from Bio.Blast import NCBIXML
from module import SerotypeHelper, JsonHelper
import logging

RESULT_FILE = "output/blast_result.xml"
DICT_FILE = "output/genome_dict.json"
BLACKLIST_FILE = "output/blacklist_genomes.json"

def create_genome_result(result_file=RESULT_FILE, dict_file=DICT_FILE, blacklist_file=BLACKLIST_FILE):
    genome_dict = defaultdict(list)
    blacklist_genome = []
    blast_result = NCBIXML.parse(open(RESULT_FILE))
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
            if genome_desc in blacklist_genome:
                continue
            genome_serotypes = SerotypeHelper.getSerotypes(genome_desc)
            is_mismatch = SerotypeHelper.isMismatch(allele_serotype, genome_serotypes)
            is_same_class = SerotypeHelper.isSameClass(allele_serotype, genome_serotypes)
            for hsp in alignment.hsps:
                percent_identity = hsp.identities / iteration.query_length
                if percent_identity < 0.97:
                    continue
                if percent_identity == 1 and is_mismatch:
                    blacklist_genome.append(
                        genome_desc + 
                        " evident: blast_result.xml {} iter {} alignment".format(iteration.query_id, alignment.hit_id)
                    )
                    continue
                if not is_same_class:
                    continue
                new_entry = {
                    "allele_desc": allele_desc,
                    "identity": percent_identity,
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
        genome_dict_sorted[genome_desc] = sorted(alignments, key=lambda alignment: alignment["identity"], reverse=True)
        
    JsonHelper.write_to_json(genome_dict_sorted,DICT_FILE)
    JsonHelper.write_to_json(blacklist_genome, BLACKLIST_FILE)

def addUsefulAllele():
    genome_dict = JsonHelper.read_from_json(DICT_FILE)
    if not genome_dict:
        logging.info("no genome_dict "+DICT_FILE+" found")
        return
    for genome_dict, alignments in genome_dict.items():
        if "wzx":
            pass
    # TODO