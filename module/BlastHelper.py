import subprocess
import os
GENOME_DIR = "data/modified_genome/"
ALLELE_FILE = "data/EcOH.fasta"
DB_FILE = "data/DB/SerotypedGenome.fasta"
RESULT_FILE = "output/blast_result.xml"
def makeBlastDB(genome_dir=GENOME_DIR, allele_file=ALLELE_FILE, db_file=DB_FILE, result_file=RESULT_FILE):
    OUTPUT_DIR = os.path.split(DB_FILE)[0]
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    OUTPUT_DIR2 = os.path.split(RESULT_FILE)[0]
    if not os.path.exists(OUTPUT_DIR2):
        os.makedirs(OUTPUT_DIR2)
    cmd = str("cat "+ GENOME_DIR +"* >> " + DB_FILE)
    cmd2 = str("makeblastdb -dbtype nucl -in " + DB_FILE)
    subprocess.call(cmd, shell=True)
    subprocess.call(cmd2, shell=True)

def blastn():
    cmd = str("blastn -query "+ALLELE_FILE+" -db "+ DB_FILE +" -out " + RESULT_FILE + " -outfmt 5")
    subprocess.call(cmd, shell=True)