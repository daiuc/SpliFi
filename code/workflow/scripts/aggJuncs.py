'''Aggregate junctions

Description: Aggregate junctions into per individual for GTEx

'''


def main(args):

    metadata, inDir, outDir, tissue = args.metadata, args.inDir, args.outDir, args.tissue

    sys.stderr.write(f'Start aggregating junctions in GTEx: {tissue} ... {datetime.now()}\n')
    sys.stderr.write(f'This script assumes all junction files {inDir} have the exact same coordinates (col1,2,3,4,6)\n')

    # GTEx metadata for all samples present in the junction file
    Gtex_Metadata = pd.read_csv(metadata, sep='\t')

    # remove white space and fix parentethies in tissue name
    fixChr = str.maketrans({" ": None, "(": "_", ")": "_"})
    Gtex_Metadata[['SMTS', 'SMTSD']] = Gtex_Metadata[['SMTS', 'SMTSD']].apply(lambda x: x.str.translate(fixChr))
    Gtex_Metadata.set_index(['SAMPID', 'SMTS', 'SMTSD', 'SUBJID'], drop=False, inplace=True)

    if tissue in list(Gtex_Metadata.SMTSD):
        sys.stderr.write(f'Use SMTSD (sub tissue category) to select tissue type: {tissue}\n')
        df = Gtex_Metadata.query('SMTSD in @tissue').reset_index(drop=True).drop_duplicates()
    else:
        sys.stderr.write(f'tissue type: {tissue} not found in SMTSD. Exiting..\n')
        exit(1)
    
    samps = {} # {subjid: [sampid, sampid, ...]}
    for k,v in df.groupby('SUBJID'):
        samps[k] = list(v.SAMPID)

    sys.stderr.write(f'aggregate junction files for Tissue type: {tissue}\n')
    sys.stderr.write(f'Total {len(samps.keys())} individuals, {len([x for l in list(samps.values()) for x in l])} samples ...\n')

    for k, v in samps.items():
        junc_files = [inDir + '/' + x + '.tsv.gz' for x in v] # junc files
        out_file = outDir + '/' + k + '.tsv'
        cts = [] # list of list for counts
        for j in junc_files:
            with gzip.open(j) as f:
                lines = [ln.decode().strip().split() for ln in f.readlines()]
                ct = [int(x[4]) for x in lines]
                cts.append(ct)
        cts = reduce(lambda x,y: [a+b for a,b in zip(x,y)], cts) # sum reads by SUBJID if having multiple SAMPID
        buf = ['\t'.join([lines[i][0], lines[i][1], lines[i][2], 
                            lines[i][3], str(cts[i]), lines[i][5]]
                        ) + '\n' for i in range(len(lines))]
        with open(out_file, 'w') as outf:
            outf.writelines(buf)
            sys.stderr.write(f'Aggregation done for {k} at {datetime.now()}\n')
    sys.stderr.write(f'Aggregation done for {tissue} at {datetime.now()}\n')


# if main then run
if __name__ == '__main__':

    from functools import reduce
    import argparse
    import sys
    import gzip
    from datetime import datetime

    import pandas as pd


    parser = argparse.ArgumentParser(description='Aggregate junctions of multiple samples into per individual for GTEx')
    parser.add_argument('--metadata', type=str, help='GTEx metadata for all samples present in the junction file')
    parser.add_argument('--inDir', type=str, help='directory that houses input junction files')
    parser.add_argument('--outDir', type=str, help='directory that houses aggregated junction files')
    parser.add_argument('--tissue', type=str, help='Tissue type')
    args = parser.parse_args()

    main(args)
