import subprocess
import os
INPUT_DIR = "data/modified_genome/"
INPUT_File = "data/EcOH.fasta"
OUTPUT_FILE = "data/DB/SerotypedGenome.fasta"
OUTPUT_FILE2 = "output/blast_result.xml"
def makeBlastDB():
    OUTPUT_DIR = os.path.split(OUTPUT_FILE)[0]
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    OUTPUT_DIR2 = os.path.split(OUTPUT_FILE2)[0]
    if not os.path.exists(OUTPUT_DIR2):
        os.makedirs(OUTPUT_DIR2)
    cmd = str("cat "+ INPUT_DIR +"* >> " + OUTPUT_FILE)
    cmd2 = str("makeblastdb -dbtype nucl -in " + OUTPUT_FILE)
    subprocess.call(cmd, shell=True)
    subprocess.call(cmd2, shell=True)

def blastn():
    cmd = str("blastn -query "+INPUT_File+" -db "+ OUTPUT_FILE +" -out " + OUTPUT_FILE2 + " -outfmt 5")
    subprocess.call(cmd, shell=True)