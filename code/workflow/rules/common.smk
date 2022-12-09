import gzip
import os
import glob
import pandas as pd



# Limit to autosomes
CHROMS = ['chr'+str(i) for i in range(1,23)]


#-------------------- Geuvadis --------------------

# Get Geuvadis metadata lookup table, key cols are Population ID, Sample ID, Run ID
Geuvadis_Metadata = pd.read_csv(config['Dataset']['Geuvadis']['Metadata'], sep='\t')
Geuvadis_Metadata.set_index(['Pop_id', 'Sample', 'run_id'], drop=False, inplace=True)

# GTEx metadata for all samples present in the junction file
Gtex_Metadata = pd.read_csv(config['Dataset']['GTEx']['Junc_meta'], sep='\t')
Gtex_Metadata.set_index(['SAMPID', 'SMTS', 'SMTSD', 'SUBJID'], drop=False, inplace=True)


# Get Geuvadis samples that have genotype in 1KGP
# For now, use only samples with 1-to-1 Sample-to-ERR
# For production, use the full linked samples
Geuvadis_Linked_SampleIDs = []
with open(config['Dataset']['Geuvadis']['Linked_1to1_SampleIDs']) as f:
    Geuvadis_Linked_SampleIDs = [s.strip() for s in f.readlines()]

def getERR(sample_id, lookup_df):
    ERR = lookup_df.query('Sample == @sample_id').run_id
    return list(set(ERR))

Sample_ERR_Dict = dict(zip(Geuvadis_Linked_SampleIDs, 
    [getERR(x, Geuvadis_Metadata) for x in Geuvadis_Linked_SampleIDs]))





