from module import DownloadMetadata, FilterGenome, RetrieveGenome, BlastHelper, SerotypeHelper, ResultAnalyzer
import os

def main():
    MORIA_DB_DIR = '/home/sam/moria/enterobasdb/'
    EXP_FILE = "data/experiment.json"
    FILTERED_FILENAMES_FILE = 'data/genome_filenames_filtered.json'
    GENOME_DIR = "data/filtered_genome/"
    ALLELE_FILE = "data/EcOH.fasta"
    BLAST_DB = "data/DB/SerotypedGenome.fasta"
    BLAST_RESULT = "output/blast_result.xml"

    if not os.path.exists(MORIA_DB_DIR):
        raise Warning("Database is not found!")
    # 1. Get the metadata from enterobase.warwick.ac.uk
    DownloadMetadata.download_metadata(EXP_FILE)
    # 2. Get a list of all genomes filename that has corresponding serotyped meta file
    FilterGenome.filter_genome(MORIA_DB_DIR, FILTERED_FILENAMES_FILE)
    # 3. copy all qualified genome data in from database
    # and Add serotype and assembly barcode to header for easier mapping later on
    RetrieveGenome.retrieve_genome(MORIA_DB_DIR, FILTERED_FILENAMES_FILE, GENOME_DIR)
    # 5. Create blast database from extracted+modified genome files
    BlastHelper.makeBlastDB(GENOME_DIR, BLAST_DB)
    # 6. Query the database with all the allele file we have
    BlastHelper.blastn(ALLELE_FILE, BLAST_DB, BLAST_RESULT)
    # 7. Create a json formatted serotype dictionary of all known/confident allele sequences
    ResultAnalyzer.create_genome_result(BLAST_RESULT)
    SerotypeHelper.initialize_dict(ALLELE_FILE)
    ResultAnalyzer.addUsefulAllele()

if __name__ == '__main__':
    main()