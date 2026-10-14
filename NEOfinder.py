import argparse
import subprocess
import re
import textwrap
import os
import sys
import uuid
import pdb
try:
    from Bio import SeqIO
except:
    print("SeqIO module is not installed! Please install SeqIO and try again.")
    sys.exit()

try:
    import tqdm
except:
    print("tqdm module is not installed! Please install tqdm and try again.")
    sys.exit()

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', required=True, type=str, dest='filename', help='Specify a fastafile.\n')
parser.add_argument('-o', '--output', required=True, dest='out', help='Specify a output fasta file name.\n')

results = parser.parse_args()
filename = results.filename
out = results.out

out_aa = out.rsplit(".",1)[0] + "_aa.fasta"
out_dna = out.rsplit(".",1)[0] + "_dna.fasta"

ori_stdout = sys.stdout

tmp = str(uuid.uuid4())
tmp_pre = str(uuid.uuid4())
os.system('> ' + out_aa)
os.system('> ' + out_dna)

stop_codon = ["TAG", "TAA", "TGA"]

proc = subprocess.Popen("grep -c '>' " + filename, shell=True, stdout=subprocess.PIPE, text=True)
length = int(proc.communicate()[0].split('\n')[0])

with tqdm.tqdm(range(length)) as pbar:
    pbar.set_description('Finding frames without stop codon...')
    for record in SeqIO.parse(filename, "fasta"):
        sys.stdout = ori_stdout
        pbar.update()
        for i in range(3):
            if "*" not in str(record.seq.reverse_complement()[i:].translate()):
                sys.stdout = ori_stdout
                f = open(tmp_pre, 'a')
                sys.stdout = f
                print(record.format("fasta"))
                break
sys.stdout = ori_stdout
print("Predicting ORFs from concatenated cDNA sequences...")
os.system(f"ORFfinder -in {tmp_pre} -out {tmp} -s 1 -strand minus -outfmt 1")

proc = subprocess.Popen("grep -c '>' " + tmp, shell=True, stdout=subprocess.PIPE, text=True)
length = int(proc.communicate()[0].split('\n')[0])

with tqdm.tqdm(range(length)) as pbar:
    pbar.set_description('Removing ORFs with stop codon...')
    for record in SeqIO.parse(tmp, "fasta"):
        pbar.update()
        ORF_no = record.description.split()[1].split("_",1)[0]
        acc = record.description.split()[1].split("_", 1)[1].split(":")[0]
        pos = "-".join(record.description.split()[1].split("_", 1)[1].split(":")[1:])
        if record.seq[-3:] not in stop_codon and len(str(record.seq)) >= 240:
            f = open(out_aa, 'a')
            sys.stdout = f
            print(f">{acc} {ORF_no} {pos}")
            print(re.sub("(.{60})", "\\1\n", str(record.seq.translate()), 0, re.DOTALL))
            sys.stdout = ori_stdout
            f = open(out_dna, 'a')
            sys.stdout = f
            print(f">{acc} {ORF_no} {pos}")
            print(re.sub("(.{60})", "\\1\n", str(record.seq), 0, re.DOTALL))
            sys.stdout = ori_stdout
os.system(f'rm {tmp} {tmp_pre}')
