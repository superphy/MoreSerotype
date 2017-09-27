import os

DB_DIR = '/home/sam/moria/enterobase_db/'
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
TEMP_DIR = os.path.join(ROOT_DIR, 'temp')
DATA_DIR = os.path.join(ROOT_DIR, 'data')
OUTPUT_DIR = os.path.join(ROOT_DIR, 'output')
SEROTYPED_ALLELE = os.path.join('data', 'EcOH.fasta')
GENOME_DIR = os.path.join(TEMP_DIR, 'genomes')
BLAST_DB = os.path.join(TEMP_DIR, 'blast_db', 'serotyped_blastdb')
GENE_LIST = ['wzx', 'wzy', 'wzm', 'wzt', 'fliC', 'fllA', 'flkA', 'flmA', 'flnA']

def make_directories():
    directories_to_make = [TEMP_DIR, DATA_DIR, OUTPUT_DIR, GENOME_DIR, BLAST_DB]
    for dir in directories_to_make:
        if not os.path.exists(dir):
            os.makedirs(dir)
