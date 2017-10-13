'''
check the +/- of the 10 E. coli specific regions for each of the downloaded enterobase genomes,
and compile that as a table.

eg. `name` `marker 1` `marker 2` `marker N`
'''
import csv
import logging
import os
import subprocess
import sys
import tempfile
from tqdm import tqdm
from multiprocessing import pool

logging.basicConfig(filename='gene_presence_analysis.log', level=logging.DEBUG)

MARKER_FILE = 'ecoli_specific_markers.fasta'

def analyze_gene_presence(fasta_files):
    '''
    Args:
        list[str]: list of paths of files to be added to gene presence report
    Returns:
        str: path of csv file that contains presence/absence of each marker on each genome
    '''
    invalid_file_count = 0
    total_hit = 0
    curr_valid_file_count = 0
    csv_file = 'gene_presence_report.csv'
    fieldnames = [
        "genome_name",
        "1436893830000|3159571",
        "1436893909000|3159808",
        "2873786891000|3159389",
        "2873787160000|3160196",
        "4310679577000|3158082",
        "4310679772000|3158667",
        "4310679831000|3158844",
        "4310680254000|3160113",
        "4310680315000|3160296",
        "4310680399000|3160548"
    ]
    with open(csv_file, 'w') as csv_handle:
        csv_writer = csv.DictWriter(csv_handle, fieldnames=fieldnames)
        csv_writer.writeheader()
    for fasta_file in tqdm(fasta_files):
        presence_output, num_hit = get_presence_output(fasta_file)
        if num_hit is -1 or num_hit is -2:
            invalid_file_count += 1
            if num_hit is -2:
                logging.warning('Blast error for %s', fasta_file)
            continue
        curr_valid_file_count += 1
        total_hit += num_hit
        msg = str("Invalid file count: %d, Average Hits: %d/%d=%.2f"
            %(invalid_file_count, total_hit, curr_valid_file_count, total_hit/curr_valid_file_count))
        tqdm.write(msg)
        logging.debug(msg)
        msg = str('%s'%list(presence_output.values()))
        tqdm.write(msg)
        logging.info(msg)
        # need lock for this part
        with open(csv_file, 'a') as csv_handle:
            csv_writer = csv.DictWriter(csv_handle, fieldnames=fieldnames)
            csv_writer.writerow(presence_output)

def get_presence_output(fasta_file):
    '''
    Args:
        fasta_file(str): fasta_file to be queried by marker database
    Returns:
        dict: dictionary of presence output
    '''
    if os.path.getsize(fasta_file) < 1000:
        logging.debug('%s is not a valid fasta file (too small)'%fasta_file)
        return None, -1
    with tempfile.TemporaryDirectory() as temp_dir:
        genome_name = os.path.splitext(os.path.basename(fasta_file))[0]
        blast_db = os.path.join(temp_dir, genome_name)
        cmd = [
            "makeblastdb",
            "-in", fasta_file,
            "-dbtype", "nucl",
            '-title', genome_name,
            '-out', blast_db
        ]
        output = run_subprocess(cmd)

        cmd2 = [
            'blastn',
            '-query', MARKER_FILE,
            '-db', blast_db,
            '-perc_identity', "90",
            '-qcov_hsp_perc', "90",
            '-max_target_seqs', '1', # this number needs to be greater than number of genome
            '-max_hsps', '1', # we only want to know hit/no hit
            '-outfmt', '6 qseqid'
        ]
        blast_result = run_subprocess(cmd2)
        if blast_result is None:
            msg = str('Blast error')
            logging.debug(msg)
            return None, -2
        new_entry = {
            "genome_name": genome_name,
            "1436893830000|3159571":0,
            "1436893909000|3159808":0,
            "2873786891000|3159389":0,
            "2873787160000|3160196":0,
            "4310679577000|3158082":0,
            "4310679772000|3158667":0,
            "4310679831000|3158844":0,
            "4310680254000|3160113":0,
            "4310680315000|3160296":0,
            "4310680399000|3160548":0
        }
        num_hit = 0
        for line in blast_result.split('\n'):
            allele_name = '|'.join(line.strip().split('|')[1:])
            if new_entry.get(allele_name) is 0:
                new_entry[allele_name] = 1
                num_hit += 1
        return new_entry, num_hit

def run_subprocess(cmd):
    output = None
    try:
        completed_process = subprocess.run(
            cmd,
            universal_newlines=True,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        logging.debug(completed_process.stderr)
        logging.info(completed_process.stdout)
        output = completed_process.stdout
    except subprocess.CalledProcessError as err:
        logging.warning("Error: %s %s", err.stderr, err.stdout)
    return output

def main(enterobase_dir = '/home/sam/moria/enterobase_db'):
    basenames = os.listdir(enterobase_dir)
    filenames = []
    for basename in basenames:
        filename = os.path.join(enterobase_dir, basename)
        filenames.append(filename)
    analyze_gene_presence(filenames)
    # analyze_gene_presence(['/home/sam/Projects/MoreSerotype/mini_experiment/ESC_HA8367AA_AS.fasta', '/home/sam/Projects/MoreSerotype/mini_experiment/ESC_OA1793AA_AS.fasta', '/home/sam/Projects/MoreSerotype/mini_experiment/Straphylococcus.fasta'])

if __name__ == '__main__':
    main()
