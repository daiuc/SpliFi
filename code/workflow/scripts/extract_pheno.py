#!/usr/bin/env python
'''
Prepare Geuvadis dataset phenotype

Goals: 
    1.  change ERR id to sample id
    2.  select columns: 
        - only samples that have a genotype
        - only samples that have 1 run
    3.  select rows: only noisy splicing (ones with "*")
    4.  output:
            - count table of selected rows and columns
            - list of sample IDs selected (used for prep genotype)
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
    help="input noise count file")

parser.add_argument("-O", "--outcount", dest="outcount",
    type=str, required=True,
    help="output count file")

parser.add_argument("--samplelist", dest="samplelist",
    type=str, required=True,
    help="output sample list file")

parser.add_argument("-L", "--lookuptable", dest="lookuptable",
    type=str, required=True,
    help="Geuvadis metadata look up table e.g. ERR <-> Sample")   

parser.add_argument("-S", "--subset", dest="subset",
    type=str, required=True,
    help="Geuvadis samples with only 1 ERR ID that has genotype")


def getERR(sample_id, lookup_df):
    ERR = lookup_df.query('Sample == @sample_id').run_id
    return list(set(ERR))

def getSample(ERR, lookup_df):
    sampleID = lookup_df.query('run_id == @ERR').Sample
    return sampleID[0]


options = parser.parse_args()

Geuvadis_Metadata = pd.read_csv(options.lookuptable, sep='\t')
Geuvadis_Metadata.set_index(['Pop_id', 'Sample', 'run_id'], 
    drop=False, inplace=True)

Geuvadis_Linked_SampleIDs = []
with open(options.subset) as f_subset:
    Geuvadis_Linked_SampleIDs = [s.strip() for s in f_subset.readlines()]


N = 0
fin = gzip.open(options.inputfile)
if os.path.exists(options.outcount):
    print(f"{options.outcount} exists. Removing it first.")
    os.remove(options.outcount)
    fout = gzip.open(options.outcount, 'wt')
    fout2 = open(options.samplelist, 'w')
else:
    fout = gzip.open(options.outcount, 'wt')
    fout2 = open(options.samplelist, 'w')

for ln in fin:
    # if N > 100:
    #     break
        
    if N == 0:
        headers = ln.decode().split()
        headers = [x.split(".")[0] if x != "chrom" else x for x in headers]
        new_headers = [getSample(c, Geuvadis_Metadata) if c != "chrom" else c for c in headers]
        select_cols = [c for c in new_headers if c in ['chrom'] + Geuvadis_Linked_SampleIDs]
        select_idx = [new_headers.index(x) for x in select_cols]
        select_samples = [c for c in select_cols if c != "chrom"]
        fout.write(' '.join(select_cols) + '\n') # write count file header
        fout2.write('\n'.join(select_samples)) # write sample ID list
        
    if N > 0:
        ln = ln.decode().split()
        select_data = [ln[idx] for idx in select_idx]
        if '*' in select_data[0]:
            fout.write(' '.join(select_data) + '\n') # write count file data

    N += 1

fin.close()
fout.close()
