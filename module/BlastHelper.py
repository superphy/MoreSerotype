import subprocess
import os
def makeBlastDB(genome_dir, db_file):
    output_dir = os.path.split(db_file)[0]
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    cmd = str("cat "+ genome_dir +"* >> " + db_file)
    cmd2 = str("makeblastdb -dbtype nucl -in " + db_file)
    subprocess.call(cmd, shell=True)
    subprocess.call(cmd2, shell=True)

def blastn(allele_file, db_file, result_file):
    result_dir = os.path.split(result_file)[0]
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
    cmd = str("blastn -query "+allele_file+" -db "+ db_file +" -out " + result_file + " -outfmt 5")
    subprocess.call(cmd, shell=True)