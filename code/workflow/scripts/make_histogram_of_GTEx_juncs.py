#!/usr/bin/env python

'''Compute histogram of GTEx junction counts

Compute histogram with bins of [0, 5, 10, 15, 20] of GTEx
junction counts. Because the junction counts file takes
over 25G of memory, I compute this histogram first. 

The purpose is to use this histogram data to verify that 
this file does indeed include lowly expressed junctions.

'''

import os
import gzip
import numpy as np
import pandas as pd


__author__    = "Chao Dai"
__email__     = "chaodai@uchicago.edu"
__status__    = "Development"
__version__   =  "v0.0.1"


if snakemake:
    juncs = snakemake.input[0]
    outfile = snakemake.output[0]
else:
    juncs = "../code/resources/GTEx/juncs/GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct.gz"

# data store a summary of each line
# first 4 numbers are histogram bins, each represent number of samples 
# falling within each bin. The 4 bins are 4 equal bins within [0, 20].
# the last 3 numbers are min, median, and max of the row. 
data = []
i = 0
with gzip.open(juncs) as f:
    for ln in f:
        if i % 5000 == 0:
            print(f"Processed {(i // 5000) * 5000} records.")
        if i == 2:
            header = ln.decode().split('\t')
        if i > 2:
            row = ln.decode().split('\t')
            keys = row[0:2]
            values = [int(x) for x in row[2:len(row)]]
            bins = np.histogram(values, bins=4, range=(0, 20))[0]
            bins = [d for d in bins]
            bins.append(np.min(values))
            bins.append(np.median(values))
            bins.append(np.max(values))
            data.append(bins)
        i += 1

# write to file
df = pd.DataFrame(np.array(data))
df.columns = ['bin_0-5', 'bin_5-10', 'bin_10-15', 'bin_15-20', 'Min', 'Med', 'Max']
df.to_csv(outfile, sep='\t', header=True, index=False)
print(f"Done. Processed {i-1} lines of records.")