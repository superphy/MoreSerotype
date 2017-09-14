from collections import defaultdict
from collections import OrderedDict
import json
import re
from Bio.Blast import NCBIXML
RESULT_FILE = "output/blast_result.xml"
OUTPUT_FILE = "output/close_mismatch_list.json"
OUTPUT_FILE2 = "output/exact_mismatch_dict.json"
OUTPUT_FILE3 = "output/exact_mismatch_list.json"

def parse_serotype_value_helper(letter, str):
    '''
    find the numeric labels following a letter from a string and return it
    return "" if nothing is found
    '''
    regex = re.compile("("+letter+"\d{1,3})(?!\d)")
    results = regex.findall(str)
    num_of_result = len(results)
    if num_of_result > 0:
        return results[0][1:]
    return ""

def parse_serotype_value(str):
    return {"H": parse_serotype_value_helper("H", str),
            "O": parse_serotype_value_helper("O", str)}

def get_serotype_classification(serotype):
    if serotype["H"]:
        return "H"
    if serotype["O"]:
        return "O"
    print("no classfication.")
    return ""

def is_serotype_mismatch(allele_serotype_classification, allele_serotype, genome_serotype):
    allele_serotype_id = allele_serotype[allele_serotype_classification]
    genome_serotype_id = genome_serotype[allele_serotype_classification]
    if allele_serotype_id:
        if genome_serotype_id:
            if allele_serotype_id != genome_serotype_id:
                return True
    else:
        print("incorrect classification")
    return False

def parse_result():
    close_mismatch_list = []
    exact_mismatch_list = []
    mismatch_count = 0
    useless_allele_count = 0
    pair_count = 0
    close_mismatch_dict = defaultdict(int)
    exact_mismatch_dict = defaultdict(int)
    # Save the top results for each query search based on percent identity

    for iteration in NCBIXML.parse(open(RESULT_FILE)):
        # Boolean value to track wheter to skip to next iteration
        is_skippable = False
        # get required attribute at this level
        allele_header = iteration.query

        looked_genome_barcode = defaultdict(bool)
        if not iteration.alignments:
            # no alignment found for this query
            continue
        iteration.alignments.sort(key=lambda align: max(hsp.identities for hsp in align.hsps),
                                  reverse=True)
        # Sort by percentage identity which is simplified to just identities
        # because query length is same for each alignment
        # eg hsp.identities / iteration.query_length == hsp.identities when sorting
        for alignment in iteration.alignments:
            pair_count += 1
            if is_skippable:
                break
            # get required attribute at this level
            genome_header = alignment.hit_def
            query_length = iteration.query_length

            # search of serotype tag with regular expression
            allele_serotype = parse_serotype_value(allele_header)
            genome_serotype = parse_serotype_value(genome_header)
            allele_serotype_classification = get_serotype_classification(allele_serotype)
            if allele_serotype_classification=="":
                # the result is meaningless if allele has no serotype
                print(allele_header)
                is_skippable = True
                useless_allele_count+=1
                continue

            # Prevent checking the same genome twice
            if looked_genome_barcode[genome_header] is True:
                continue
            looked_genome_barcode[genome_header] = True

            for hsp in alignment.hsps:
                percent_identity = hsp.identities / query_length

                # if percentage identity is too low, we skip to next query
                if percent_identity < .97:
                    is_skippable = True
                    break
                if is_serotype_mismatch(allele_serotype_classification, allele_serotype, genome_serotype):
                    new_result = {
                        "Allele Header": allele_header,
                        "Genome Header": genome_header,
                        "Query Length": query_length,
                        "Alignment Length": hsp.align_length,
                        "Allele Serotype": allele_serotype,
                        "Genome Serotype": genome_serotype,
                        "Score": hsp.score,
                        "Percent_identity": str(percent_identity * 100) + "%",
                        "Align Seq": {
                            "query": hsp.query,
                            "match": hsp.match,
                            "sbjct": hsp.sbjct
                        }
                    }
                    if percent_identity == 1:
                        mismatch_count+=1
                        exact_mismatch_list.append(new_result)
                        exact_mismatch_dict[allele_serotype_classification+allele_serotype[allele_serotype_classification]]+=1
                    else:
                        close_mismatch_list.append(new_result)
                        close_mismatch_dict[allele_serotype_classification+allele_serotype[allele_serotype_classification]]+=1
    with open(OUTPUT_FILE, 'w') as outfile:
        json.dump(close_mismatch_list, outfile, indent=4, separators=(',', ': '))
    
    # store mismatch summary as a sorted dictionary
    exact_mismatch_dict_sorted = OrderedDict(sorted(exact_mismatch_dict.items(), key=lambda t: t[1]))
    with open(OUTPUT_FILE2, 'w') as outfile:
        json.dump(exact_mismatch_dict_sorted, outfile, indent=4, separators=(',', ': '))
    with open(OUTPUT_FILE3, 'w') as outfile:
        json.dump(exact_mismatch_list, outfile, indent=4, separators=(',', ': '))
    print("Operation completed.")
    print("{} mismatches encountered".format(mismatch_count))
    print("{} useless alleles encountered".format(useless_allele_count))
    print("{} potential candidate identified and stored in close_mismatch_list.json".format(len(close_mismatch_list)))
    print("{} total alignments".format(pair_count))