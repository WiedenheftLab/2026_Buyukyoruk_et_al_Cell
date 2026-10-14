#!/usr/bin/python

import os
import subprocess

proc = subprocess.Popen("pwd", shell=True, stdout=subprocess.PIPE, )
dir = (proc.communicate()[0].split('\n')[0])

filenames = os.listdir(dir)
for filename in filenames:
	if '.fasta' in filename:
		os.system('prodigal -i ' + filename + " -a " + filename.split(".fasta")[0] + "_prodigal_aa.fasta")

