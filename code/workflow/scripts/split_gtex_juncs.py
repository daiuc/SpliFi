#!/usr/bin/env python
'''
Munge GTEX junction file

Goals: 
    1.  split munged GTEX junction file columns to each individual file
        with the following format:
        chr1    12058   12178   .       0       + 
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

parser.add_argument("-O", "--outputfile", dest="outputfile",
    type=str, required=True,
    help="output count file")

parser.add_argument("-S", "--sampleid", dest="sampleid",
    type=str, required=True,
    help="GTEx SAMID to select")

options = parser.parse_args()

samid = options.sampleid # sample id to select

fin = gzip.open(options.inputfile)
fout = gzip.open(options.outputfile, 'wt')

print('split munged output..')


N = 0
for ln in fin:
    # if N > 10000: # diagnosis
    #     break
    
    if N > 100 and N % 5000 == 0:
        print(f"Written {N} lines..")

    if N == 0: # read header row
        header = ln.decode().split()

        try:
            idx = header.index(samid)
        except ValueError:
            print("Sample ID is not in the file! Exit.")
            exit(1)

    if N > 0: # data rows
        data = ln.decode().split()
        chrom, start, end, dot, strand = data[0:5]
        count = data[idx]
        data = [chrom, start, end, dot, count, strand]

        # write data
        fout.write('\t'.join(data) + '\n')
      
    N += 1

print(f"Written {N} lines. Done")
fin.close()
fout.close()


                




















