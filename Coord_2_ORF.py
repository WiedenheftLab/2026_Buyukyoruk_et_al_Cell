#!/software/mambaforge/envs/Murat_scripts/bin/python

import argparse
import os
import sys
import subprocess
import pandas as pd

try:
    import tqdm
except ImportError, e:
    print "tqdm module is not installed! Please install tqdm and try again."
    sys.exit()

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True, type=str, dest='filename',
                    help='Specify a fastafile to fetch regions from.\n')
parser.add_argument('-o', '--output', required=True, dest='out',
                    help='Specify a output file name.\n')
parser.add_argument('-d', '--data', required=True, dest='data',
                    help='Specify a dataframe with accession, start, stop, strand info in that order.\n')

results = parser.parse_args()
filename = results.filename
out = results.out
data = results.data

df = pd.read_csv(filename, sep='\t', lineterminator='\n')
df_extract = pd.read_csv(data, sep='\t', lineterminator='\n')
final_res = pd.DataFrame()

for index,row in df_extract.iterrows():
    acc = row['Accession']
    start = row['start']
    end = row['stop']
    strand = row['strand']
    if strand == 1:
        df_found = df[df['Accession'].isin([acc]) & df['strand'].isin([strand]) & df['end'].isin([end])]
    else:
        df_found = df[df['Accession'].isin([acc]) & df['strand'].isin([strand]) & df['start'].isin([start])]
    if not df_found.empty:
        final_res = pd.concat([final_res, df_found['ORF']]).reset_index(drop=True)

final_res.to_csv(out, sep='\t', index=False, header=False)




