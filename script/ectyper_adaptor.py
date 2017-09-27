import os
import sys
from collections import defaultdict

import Bio
import Bio.Seq
import Bio.SeqIO
import Bio.SeqRecord

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

import definitions
import module.JsonHelper


def adapt_ectyper(serotype_dict_f):
    serotype_dict = \
        module.JsonHelper.read_from_json(serotype_dict_f)
    ectyper_dict = defaultdict(dict)
    ectyper_data = []
    for serotype, allele_list in serotype_dict.items():
        for allele in allele_list:
            # information for ectyper_dict entry
            serotype_class = serotype[0]
            allele_id = '-'.join([serotype, str(allele['num'])])
            gene = allele['gene']
            serotype = serotype

            ectyper_dict_entry = {
                "gene" : gene,
                "allele" : serotype
            }
            ectyper_dict[serotype_class][allele_id] = ectyper_dict_entry
            # information for ectyper_data entry
            allele_id = allele_id
            allele_seq = allele['seq']

            ectyper_data_entry = Bio.SeqRecord.SeqRecord(
                Bio.Seq.Seq(
                    allele['seq'],
                    Bio.Alphabet.IUPAC.IUPACUnambiguousDNA),
                id = allele_id,
                description = "")
            ectyper_data.append(ectyper_data_entry)
    # store output
    module.JsonHelper.write_to_json(
        ectyper_dict,
        os.path.join(definitions.OUTPUT_DIR, 'ectyper_dict.json')
    )
    Bio.SeqIO.write(
        ectyper_data,
        os.path.join(
            definitions.OUTPUT_DIR,
            'ectyper_data.fasta'
        ),
        'fasta'
    )
    allele_count = len(ectyper_data)
    print("%d entries compatible to ectyper now" %allele_count)

def main():
    adapt_ectyper(os.path.join(definitions.OUTPUT_DIR, 'serotype_dict.json'))

if __name__ == '__main__':
    main()
