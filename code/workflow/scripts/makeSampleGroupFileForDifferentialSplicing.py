'''
This script is used to make a sample group file for differential splicing analysis.
'''


import sys
import os
import gzip
import argparse

def writeSampleGroupFile(fname, outname):
    sys.stdout.write(f'Writing sample group file: "{outname}" from "{fname}"\n')
    with gzip.open(fname, 'rt') as f:
        header = f.readline().split()
        cols = header[1:]
        groups = [f.split('.')[0] for f in cols]
        sys.stdout.write(f'Detected {len(set(cols))} columns(samples) from the following {len(set(groups))} groups:\n')
        for g in set(groups):
            sys.stdout.write(f'{g}: {sum([g in col for col in cols])} samples\n')
        outLines = [a + ' ' + b for a, b in zip(cols, groups)]
        with open(outname, 'w') as f:
            sys.stdout.write(f'Writing to {outname}\n')
            f.write('\n'.join(outLines))
    sys.stdout.write('Done.\n')


def mungeNumersFile(fname, outname):
    sys.stdout.write(f'Reformatting {fname} to be compatible with leafcutter_ds.R\n')
    with gzip.open(fname, 'rt') as f:
        with gzip.open(outname, 'wt') as out:
            for ln in f:
                if ln.startswith('chrom'): # header
                    header = ln.removeprefix('chrom ')
                    out.write(header)
                else:
                    rec = ln.split()
                    fixed_coord = ':'.join(rec[0].split(':')[:4])
                    rec[0] = fixed_coord
                    out.write(' '.join(rec) + '\n')
    sys.stdout.write(f'Done. Wrote to {outname}\n')


def main(args):
    infile = args.inputFile
    outfile = args.outputFile
    sampfile = args.sampGroupFile
    writeSampleGroupFile(infile, sampfile)
    mungeNumersFile(infile, outfile)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Make a sample group file for differential splicing analysis.')
    parser.add_argument('-i', dest='inputFile', help='output of leafcutter2, e.g. ds_perind.constcounts.noise_by_intron.gz')
    parser.add_argument('-o', dest='outputFile', help='Output numerators file name, e.g. ds_perind.constcounts.noise_by_intron.lf1.gz')
    parser.add_argument('-s', dest='sampGroupFile', help='sample group file')
    args = parser.parse_args()
    main(args)