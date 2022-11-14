#!/usr/bin/env python
'''
Munge GTEX junction file

Goals: 
    1.  Organize junction data such that the first 5 cols are:
        chrom, A, B, dot, counts, strand
        col6 and beyond are the counts for each sample
'''

import os
import gzip
import argparse


__author__    = "Chao Dai"
__email__     = "chaodai@uchicago.edu"
__status__    = "Development"
__version__   =  "v0.0.1"


import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-I", "--inputfile", dest="inputfile",
    type=str, required=True,
    help="input count file")

parser.add_argument("-O", "--outcount", dest="outcount",
    type=str, required=True,
    help="output count file")

parser.add_argument("-G", "--gencode", dest="gencode",
    type=str, required=True,
    help="Gencode annotation, a csv dataframe file")   

options = parser.parse_args()

# dictionary of {'gene_id': 'strand'}
anno = pd.read_csv(options.gencode)
strand_lk = dict(anno[['gene_id', 'strand']].drop_duplicates().values)
print(f"Gene ID to strand lookup table, {len(strand_lk)} records.")

fin = gzip.open(options.inputfile)
fout = gzip.open(options.outcount, 'wt')


print('Writing munged output..')
N = 0
for ln in fin:
    # if N > 15: # diagnosis
    #     break
    
    if N > 100 and N % 5000 == 0:
        print(f"Written {N} lines..")

    if N == 2: # header row
        header = ln.decode().split()
        sample_keys = header[2:len(header)]
        new_header = ['chrom', 'start', 'end', 'dot', 'strand'] + sample_keys

        # write header
        fout.write(' '.join(new_header) + '\n')

    if N > 2: # data rows
        data = ln.decode().split()
        chrom, A, B = data[0].split('_')
        try:
            strand = strand_lk[data[1]]
        except KeyError:
            strand = "*"
        new_data = [chrom, A, B, '.', strand] + data[2:len(data)]

        # write data
        fout.write(' '.join(new_data) + '\n')
    
    N += 1

print(f"Written {N} lines. Done")
fin.close()
fout.close()


                




















