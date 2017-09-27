from module import JsonHelper
import os
import logging
import definitions

log = logging.getLogger(__name__)

def filter_genome(strains_file):
    '''
    create a file which contains a list of selected genome files
    from DB that have matching metadata in strains.json.
    Return the created file.
    '''
    log.info("Start filtering genome")
    strains_data = JsonHelper.read_from_json(strains_file)
    selected_genomes = []
    genome_files = os.listdir(definitions.DB_DIR)
    genome_files_len = len(genome_files)
    for index, file in enumerate(genome_files):
        log.debug("{}/{} genome files filtered".format(index, genome_files_len))
        for strain in strains_data:
            if strain['assembly_barcode'] in file:
                selected_genomes.append(file)
    log.info("%d genome files selected.", len(selected_genomes))
    selected_genomes_file = os.path.join(definitions.TEMP_DIR, 'selected_genomes.json')
    JsonHelper.write_to_json(selected_genomes, selected_genomes_file)
    return selected_genomes_file