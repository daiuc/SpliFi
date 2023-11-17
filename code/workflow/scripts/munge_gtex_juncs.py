#!/usr/bin/env python
'''
Munge GTEX junction file

Goals: 
    1.  Organize junction data such that the first 5 cols are:
        chrom, A, B, dot, counts, strand
        col6 and beyond are the counts for each sample
    2.  strand is inferred from the Gencode annotation file. Choose gencode_v26
        because it is the same version used in the junciton file.
'''

import sys
import gzip
import argparse
from datetime import datetime


__author__    = "Chao Dai"
__email__     = "chaodai@uchicago.edu"
__status__    = "Development"
__version__   =  "v0.0.1"


import pandas as pd


parser = argparse.ArgumentParser(description="Munge GTEX junction file and infer strand based on gene ID from Gencode annotation file")
parser.add_argument("-I", "--inputfile", dest="inputfile", type=str, required=True, help="input count file")
parser.add_argument("-O", "--outcount", dest="outcount", type=str, required=True, help="output count file")
parser.add_argument("-G", "--gencode", dest="gencode", type=str, required=True, help="Gencode annotation, a csv dataframe file")   
options = parser.parse_args()

inputfile, gencode, outcount = options.inputfile, options.gencode, options.outcount

sys.stderr.write(f'## Started at {datetime.now().strftime("%Y-%m-%d %H:%M:%S")} \n')
sys.stderr.write(f'## Munge GTEx junction file: {inputfile} and use {gencode} file to infer strand by gene_id. \n')


anno = pd.read_csv(gencode, sep="\t", names=['chrom', 'A', 'B', 'gene_id', 'dot', 'strand'])
strand_lk = dict(anno[['gene_id', 'strand']].drop_duplicates().values)
sys.stderr.write(f"Gencode {gencode} Gene ID to strand lookup table, {len(strand_lk)} records.\n")


fin = gzip.open(inputfile)
fout = gzip.open(outcount, 'wt')

sys.stderr.write(f'## Write munged output to {outcount} \n')
N = 0
for ln in fin:
    if N > 100 and N % 5000 == 0:
        sys.stderr.write(f"Wrote {N} lines..\n")

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

sys.stderr.write(f'## Finished at {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}. Wrote {N} lines.\n')

fin.close()
fout.close()


                




















