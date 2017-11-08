import os
import sys
import pandas as pd
import tempfile
from multiprocessing import Manager, Process
import multiprocessing

ectyper_dir = '/home/sam/Projects/ecoli_serotyping'
sys.path.append(ectyper_dir)

from Bio import SeqIO
from tqdm import tqdm

from ectyper import speciesIdentification
genome_dir = '/home/sam/moria/enterobase_db'

report = pd.read_csv('new_report.csv', index_col='genome name')
bad_predictions = report[(report['hit count']<3) &
       (report['enterobase species'].str.find('Escherichia')!=-1) &
       report['enterobase species'].notnull()]

def worker(barcode, species_results):
    print(multiprocessing.current_process())
    genome_file = os.path.join(genome_dir, barcode+'.fasta')
    if os.path.getsize(genome_file) < 1000:
        print("%s is too small.", genome_file)
        return
    with tempfile.TemporaryDirectory() as temp_dir:
        new_file = os.path.join(temp_dir, barcode+'.fasta')
        new_fh = open(new_file,'w')
        header = '> %s\n' %barcode
        new_fh.write(header)
        for record in SeqIO.parse(genome_file, 'fasta'):
            new_fh.write(str(record.seq))
            new_fh.write('nnnnnnnnnnnnnnnnnnnn')
        new_fh.close()
        try:
            species = speciesIdentification.get_species(new_file)
        except:
            species = ''
        new_entry = {
            'barcode':barcode,
            'species':species
        }
        species_results.append(new_entry)
    return
# species_results = []
# with Manager() as manager:
#     jobs = []
#     for barcode in tqdm(bad_predictions.index):
#         worker(barcode, species_results)
#         tqdm.write(str(species_results))
def main():
    pool = multiprocessing.Pool(multiprocessing.cpu_count()//2)
    with Manager() as manager:
        L = manager.list()  # <-- can be shared between processes.
        processes = []
        for barcode in bad_predictions.index:
            pool.apply_async(worker, args=(barcode, L))
        pool.close()
        pool.join()
        pd.DataFrame(list(L)).to_csv('mash_for_lowhitecoli.csv')
        print('completed')

if __name__ == '__main__':
    main()