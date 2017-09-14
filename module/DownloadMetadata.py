import requests
from module import JsonHelper

OUTPUT_FILE = "data/strains.json"
OUTPUT_FILE2 = "data/experiment.json"

def download_metadata():
    # Downloads all the E.coli genomes from Enterobase.
    options = {
        'experiment':'assembly_stats',
        'database':'ecoli',
        'strain_query_type':'query',
        'strain_query':"serotype LIKE  '%O%' OR serotype LIKE  '%H%' OR serotype LIKE  '%K%'",
        'experiment_query':"status LIKE  '%Assembled%'",
        'operand':'AND'
    }
    r = requests.post('http://enterobase.warwick.ac.uk/get_data_for_experiment', data=options)
    d = r.json()
    strains = d['strains']
    experiment = d['experiment']
    # Retrive assembly barcode from experiment
    strains_len = len(strains)
    experiment_len = len(experiment)
    print("{} strains and {} experiment objects found".format(strains_len, experiment_len))
    print("Begin mapping assembly barcode")
    for index, row in enumerate(strains):
        print("Finding assembly barcode for {}/{} strains".format(index+1, strains_len))
        assembled = row['assembly_status']
        assembly_barcode=''
        if assembled == 'Assembled':
            id = row['id']
            for row2 in experiment:
                if row2['id'] == id:
                    assembly_barcode = row2['barcode']
            if assembly_barcode == '':
                print(' Warning: no assembly barcode found')
        strains[index]['assembly_barcode'] = assembly_barcode
    
    JsonHelper.write_to_json(strains, OUTPUT_FILE)
    JsonHelper.write_to_json(experiment, OUTPUT_FILE2)