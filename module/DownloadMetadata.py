import logging
import os
import random

import requests

import definitions
from module import JsonHelper

log = logging.getLogger(__name__)

def download_metadata():
    '''
    Downloads all the E.coli genomes from Enterobase.
    Return path to 'strains.json'
    '''
    options = {
        'experiment':'assembly_stats',
        'database':'ecoli',
        'strain_query_type':'query',
        'strain_query':"serotype LIKE  '%O%' OR serotype LIKE  '%H%' OR serotype LIKE  '%K%'",
        'experiment_query':"status LIKE  '%Assembled%'",
        'operand':'AND'
    }
    r = requests.post('http://enterobase.warwick.ac.uk/get_data_for_experiment', data=options)
    log.info("Start requesting metadata from enterobase.warwick.ac.uk.")
    d = r.json()
    strains = d['strains']
    strains = strains
    experiment = d['experiment']
    # Retrive assembly barcode from experiment
    strains_len = len(strains)
    experiment_len = len(experiment)
    log.info("{} strains and {} experiment objects found".format(strains_len, experiment_len))
    log.info("Begin mapping assembly barcode from experiment.json to strains.json")
    for index, row in enumerate(strains):
        log.debug("Finding assembly barcode for {}/{} strains".format(index+1, strains_len))
        assembled = row['assembly_status']
        assembly_barcode=''
        if assembled == 'Assembled':
            id = row['id']
            for row2 in experiment:
                if row2['id'] == id:
                    assembly_barcode = row2['barcode']
            if assembly_barcode == '':
                log.warning('No assembly barcode found for row %i', row)
        strains[index]['assembly_barcode'] = assembly_barcode
    log.info("Remapping complete. Storing file.")
    strains_file = os.path.join(definitions.TEMP_DIR, 'strains.json')
    experiments_file = os.path.join(definitions.TEMP_DIR, 'experiments.json')
    JsonHelper.write_to_json(strains, strains_file)
    JsonHelper.write_to_json(experiment, experiments_file)
    return strains_file