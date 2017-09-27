import logging
import os

import definitions
from module import (BlastHelper, DownloadMetadata, FilterGenome, JsonHelper,
                    ResultAnalyzer, RetrieveGenome, SerotypeHelper, LoggingHelper)

log = logging.getLogger(__name__)

def main():

    LoggingHelper.setup_logging()
    definitions.make_directories()
    if not os.path.exists(definitions.DB_DIR):
        log.fatal("Database is not found!", exc_info=1)
    # 1. Get the metadata from enterobase.warwick.ac.uk
    strains_file = DownloadMetadata.download_metadata()
    # 2. Get a list of all genomes filename that has corresponding serotyped meta file
    selected_genomes_file = FilterGenome.filter_genome(strains_file)
    # 3. copy all qualified genome data in from database
    # and Add serotype and assembly barcode to header for easier mapping later on
    RetrieveGenome.retrieve_genomes(strains_file, selected_genomes_file)
    # 5. Create blast database from extracted+modified genome files
    BlastHelper.makeBlastDB()
    # 6. Query the database with all the allele file we have
    blast_output = BlastHelper.blastn(definitions.SEROTYPED_ALLELE, definitions.BLAST_DB)
    # 7. Create a json formatted serotype dictionary of all known/confident allele sequences
    genome_dict_file = ResultAnalyzer.create_genome_result(blast_output)
    serotype_dict_file = ResultAnalyzer.initialize_serotype_dict()
    ResultAnalyzer.expand_serotype_dict(serotype_dict_file, genome_dict_file)
    log.info("Program completed")

if __name__ == '__main__':
    main()
