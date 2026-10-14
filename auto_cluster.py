import os
import subprocess

workdir = os.path.abspath(os.getcwd())

filenames = []

for file in os.listdir(workdir):
    if file.endswith(".fasta"):
        filenames.append(file)

filenames.sort()

uclust = "/software/usearch/usearch11.0.667_i86linux32"
os.system("mkdir originals")

for i in range(len(filenames)):
    os.system(f"{uclust} -cluster_fast " + filenames[i] + " -id 1 -centroids temp.fasta")

    proc = subprocess.Popen("grep -c '>' temp.fasta" , shell=True, stdout=subprocess.PIPE, text=True)
    length = int(proc.communicate()[0].split('\n')[0])

    if length >= 10:
        os.system(f"{uclust} -cluster_fast " + filenames[i] + " -id 0.6 -centroids " + filenames[i].split('.fasta')[0] + '_uclust_60.fasta -uc ' + filenames[i].split('.fasta')[0] + '_ucluster_60.uc')
    else:
        os.system(f"{uclust} -cluster_fast " + filenames[i] + " -id 0.8 -centroids " + filenames[i].split('.fasta')[0] + '_uclust_80.fasta -uc ' + filenames[i].split('.fasta')[0] + '_ucluster_80.uc')

    os.system('mv ' + filenames[i] + ' originals')
    os.system('mv *.uc originals')
