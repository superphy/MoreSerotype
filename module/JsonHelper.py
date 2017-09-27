import json
import os
import logging

log = logging.getLogger(__name__)

def write_to_json(data, output_filepath):
    '''
    write the data object into a json file
    '''
    output_dir = os.path.split(output_filepath)[0]
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    with open(output_filepath, 'w') as handler:
        json.dump(data, handler, indent=4, separators=(',', ': '))
        handler.close()
    log.debug("data with {} elements written to {}".format(len(data), output_filepath))

def read_from_json(input_filepath):
    '''
    read data from a json file and return it
    '''
    with open(input_filepath) as handler:
        data = json.load(handler)
        handler.close()
    log.debug("Read {} elements from {}".format(len(data), input_filepath))
    return data