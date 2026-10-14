import os

workdir = os.path.abspath(os.getcwd())

filenames = []

for file in os.listdir(workdir):
    if file.endswith(".fasta"):
        filenames.append(file)

filenames.sort()

for file in filenames:
    dir_out = file.split(".fasta")[0] + "_meme"
    os.system(f"meme {file} -dna -oc {dir_out} -mod anr -nmotifs 6 -minw 6 -maxw 40 -objfun classic -markov_order 1")


