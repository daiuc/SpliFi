#!/bin/env python



def getSamples(file: str, tissue: str, sep: str='\t') -> list: 
    """
    Get list of sample names based tissue name.

    - file: str: path to the sample metadata file, by default is tsv file.
    - tissue: str: tissue name, must be how it is named in the GTEx sample metadata file.
    - sep: str: separator used in the file, by default is ','.
    
    Returns: list: list of sample names for the specified tissue.
    """

    # remove white space and replace parenthesis with underscore
    fixChr = str.maketrans({" ": None, "(": "_", ")": "_"})

    # GTEx metadata for all samples present in the junction file
    Gtex_Metadata = pd.read_csv(file, sep=sep)
    Gtex_Metadata[['SMTS', 'SMTSD']] = Gtex_Metadata[['SMTS', 'SMTSD']].apply(lambda x: x.str.translate(fixChr))
    Gtex_Metadata.set_index(['SAMPID', 'SMTS', 'SMTSD', 'SUBJID'], drop=False, inplace=True)

    samples = Gtex_Metadata.query('SMTSD == @tissue').reset_index(drop=True)['SAMPID'].tolist()

    return samples

if __name__ == "__main__":
    
    import sys
    import gzip
    import argparse
    from datetime import datetime
    import pandas as pd

    parser = argparse.ArgumentParser(description='Extract gene expression from GTEx')
    parser.add_argument('-I', dest='input', type=str, help='GTEx gene expression file downloaded from GTEx portal')
    parser.add_argument('-M', dest='meta', type=str, help='GTEx sample metadata file, e.g. sampid-smts-smtsd-subjid.tsv')
    parser.add_argument('-O', dest='output', type=str, help='Extracted gene expression for the specified tissue')
    parser.add_argument('-T', dest='tissue', type=str, help='Tissue name, must be how it is named in the GTEx sample metadata file')

    args = parser.parse_args()

    inFile, metaFile, outFile, tissue = args.input, args.meta, args.output, args.tissue
    sys.stdout.write(f'Started script at {datetime.now()} ...\n')
    sys.stdout.write(f'Extracting gene expression for {tissue} from {inFile} to {outFile} ...\n') 

    # get list of samples for the specified tissue
    sample_names = getSamples(metaFile, tissue)
    sys.stdout.write(f'Number of samples for {tissue} is {len(sample_names)} ...\n')

    # extract gene expression for the specified tissue
    outFile = gzip.open(outFile, 'wt')
    with gzip.open(inFile, 'rt') as f:
        for i, line in enumerate(f):
            if len(line) < 20: continue
            if i == 2 or line.startswith('Name'):
                header = line.strip().split('\t')
                sample_idx = [header.index(x) for x in sample_names]
                outline = '\t'.join(header[:2] + [header[x] for x in sample_idx]) + '\n'
                outFile.write(outline)
            elif i > 2:
                line = line.strip().split('\t')
                ge = [line[x] for x in sample_idx] # gene expression
                outline = '\t'.join(line[:2] + ge) + '\n'
                outFile.write(outline)

    outFile.close()
    sys.stdout.write(f'Done! {datetime.now()}\n')



