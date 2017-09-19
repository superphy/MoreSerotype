import os
import multiprocessing
import time
import json
from threading import Thread
from queue import Queue
import requests

'''
Multi-processing version of enterobase_db. Expected runtime = enterobase_db.py runtime / # CPU core
54116 files downloaded in 66205 seconds
'''

def worker():
    while True:
        (identifier, barcode) = TASK_QUEUE.get()
        get(identifier, barcode)
        TASK_QUEUE.task_done()

INPUT_FILE = "strains.json"
OUTPUT_DIR = 'enterobase_db'
NUM_WORKERS = multiprocessing.cpu_count()
TASK_QUEUE = Queue()
THREADS = [Thread(target=worker) for _ in range(NUM_WORKERS)]
START_TIME = time.time()



def get(identifier, barcode):
    task_start_time = time.time()
    fn = OUTPUT_DIR + '/' + str(barcode) + '.fasta'
    print("Start downloading {}".format(fn))
    try:
        f = requests.get('http://enterobase.warwick.ac.uk/upload/download?assembly_id=' + str(identifier) + '&database=ecoli')
    except:
        time.sleep(10)
    if not f.text == 'No Assembly or strain record found':
        with open(fn, 'w') as fl:
            fl.write(f.text)
            end_time = time.time()
            print('{} downloaded in {:.2f} seconds'.format(fn, end_time - task_start_time))
            print('Program runtime: {:.2f} seconds'.format(end_time - START_TIME))


def getAll(strains):
    for row in strains:
        identifier = row['best_assembly']
        barcode = row['assembly_barcode']
        assembled = row['assembly_status']
        if assembled == 'Assembled':
            TASK_QUEUE.put((identifier, barcode))
    for thread in THREADS:
        thread.start()

    TASK_QUEUE.join()
    program_end_time = time.time()
    print("Program completed in {:.2f} seconds".format(program_end_time - START_TIME))


def enterobase():
    '''
    Downloads all the E.coli genomes from Enterobase.
    '''
    strains = None
    if os.path.isfile(INPUT_FILE):
        print("use existing strains file")
        with open(INPUT_FILE) as handler:
            strains = json.load(handler)
            handler.close()
    else:
        options = {
            'no_legacy':'true',
            'experiment':'assembly_stats',
            'database':'ecoli',
            'strain_query_type':'query',
            'strain_query':'all'
        }
        r = requests.post('http://enterobase.warwick.ac.uk/get_data_for_experiment', data=options)
        d = r.json()
        strains = d['strains']
        experiment = d['experiment']
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
        print("Finish mapping assembly barcode")
        with open(INPUT_FILE, 'w') as outfile:
            json.dump(strains, outfile, indent=4, separators=(',', ': '))
        print("Modified strains file stored as " + INPUT_FILE)
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    getAll(strains)


if __name__ == '__main__':
    enterobase()
