import subprocess
import os
import logging
import definitions
import tempfile

log = logging.getLogger(__name__)

def makeBlastDB():
    tempdir = tempfile.gettempdir()
    db_name = 'serotyped_blastdb'
    concatenated_genome_file = os.path.join(tempdir, db_name+".fasta")
    log.info("Concatenate all genome files for makeblastdb")
    filenames = [
        os.path.join(definitions.GENOME_DIR, filename) 
            for filename in os.listdir(definitions.GENOME_DIR)
    ]
    with open(concatenated_genome_file, 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for _ in infile:
                    outfile.write(_)
    log.info("Concatenation complete")
    makeblastdb_cmd = [
        "makeblastdb",
        "-in", concatenated_genome_file,
        "-dbtype", "nucl",
        '-title', db_name,
        '-out', definitions.BLAST_DB
    ]
    log.info("Executing command: \n%s", ' '.join(makeblastdb_cmd))
    completed_process = subprocess.run(
        makeblastdb_cmd,
        check=True,
        universal_newlines=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    log.debug("Output from makeblastdb:")
    log.debug(completed_process.stdout)
    log.debug(completed_process.stderr)
    return definitions.BLAST_DB

def blastn(query_file, db_file):
    log.info('Performing blastn search with %s on %s', definitions.SEROTYPED_ALLELE, definitions.BLAST_DB)
    output_file = db_file+'.xml'
    cmd = [
        'blastn',
        '-query', query_file,
        '-db', db_file,
        '-perc_identity', "97",
        '-qcov_hsp_perc', "90",
        '-max_target_seqs', '6000', # this number needs to be greater than number of genome
        '-max_hsps', '1',
        '-out', output_file,
        '-outfmt', '5'
    ]
    completed_process = subprocess.run(
        cmd,
        check=True,
        universal_newlines=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    if completed_process.returncode == 0:
        log.info("blastn successfully completed.")
        log.debug("Output from blastn:")
        log.debug("See %s", output_file)
        log.debug(completed_process.stderr)
        return output_file
    else:
        log.fatal("blastn did not run successfully.\n%s",
                  completed_process.stderr)
        exit(1)