from module import DownloadMetadata, FilterGenome, RetrieveGenome, BlastHelper, HeaderModifier, SerotypeHelper, ResultAnalyzer

def main():
    MORIA_DB_DIR = '/home/sam/moria/enterobase_db'
    EXP_FILE = "data/experiment.json"
    FILTERED_FILENAMES_FILE = 'data/genome_filenames_filtered.json'
    GENOME_DIR = "data/filtered_genome/"
    ALLELE_FILE = "data/EcOH.fasta"
    BLAST_RESULT = "output/blast_result.xml"

    '''
    # 1. Get the metadata from enterobase.warwick.ac.uk
    DownloadMetadata.download_metadata(EXP_FILE)
    # 2. Get a list of all genomes filename that has corresponding serotyped meta file
    FilterGenome.filter_genome(MORIA_DB_DIR, FILTERED_FILENAMES_FILE)
    # 3. copy all qualified genome data in from database
    RetrieveGenome.retrieve_genome(MORIA_DB_DIR, FILTERED_FILENAMES_FILE, GENOME_DIR)
    # 4. Add serotype and assembly barcode to header for easier mapping later on
    HeaderModifier.modify_header(GENOME_DIR, "data/modified_genome")
    # 5. Create blast database from extracted+modified genome files
    BlastHelper.makeBlastDB("data/modified_genome/", "data/DB/SerotypedGenome.fasta")
    # 6. Query the database with all the allele file we have
    BlastHelper.blastn(ALLELE_FILE, "data/DB/SerotypedGenome.fasta", BLAST_RESULT)
    # 7. Create a json formatted serotype dictionary of all known/confident allele sequences
    '''
    SerotypeHelper.initialize_dict(ALLELE_FILE, "output/serotype_dict.json")
    ResultAnalyzer.create_genome_result(BLAST_RESULT)
    ResultAnalyzer.addUsefulAllele()

if __name__ == '__main__':
    main()