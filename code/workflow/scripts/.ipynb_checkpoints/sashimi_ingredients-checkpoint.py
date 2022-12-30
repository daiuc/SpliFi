#!/usr/bin/env python

'''Prepare inputs for sashimi plots

Two main ingredients:
    - averaged bigwig files by allele type
    - linkage files 
'''

__author__    = "Chao Dai"
__email__     = "chaodai@uchicago.edu"
__status__    = "Development"
__version__   =  "v0.0.1"



import os
import numpy as np
import pysam as ps
import pyBigWig as pw


def getIndex(l, m):
    '''returns index of all occurances of v in l.
    l : list to query from
    m : value match
    '''
    return [k for k,v in enumerate(l) if v == m]


def getSamplesByAllele(bcf, snp):
    '''Get sample names by allele
    bcf : pysam VariantFile object
    snp : tuple. (chrom, position) where position is 1 based inclusive.
    
    return : dict
        key   : 0, 1, or 2, representing allele
        value : sample names that belong to each allele
    '''
    chrom, pos = snp
    pos = int(pos)
    samples = list(bcf.header.samples)
    for rec in bcf.fetch(chrom, pos - 1, pos + 1):
        if rec.pos == pos:
            gt = [x['GT'] for x in rec.samples.values()]
            gt = [a+b for a,b in gt]
            gt0 = getIndex(gt, 0)
            gt1 = getIndex(gt, 1)
            gt2 = getIndex(gt, 2)
    
    return {0:[samples[x] for x in gt0],
            1:[samples[x] for x in gt1],
            2:[samples[x] for x in gt2]}


def refineSamples(d, l, bwfiles):
    '''Ensure samples appear in both genotype and phenotypes
    
    d       : dict. Dictionary of genotype and sample names, as in output of
              getSamplesByGenotype.
    l       : list. A sample name list. Sample names in genotype dict must be also 
              in this list.
    bwfiles : list of bigwig files (path).
    
    return  : 
        - gt, refined sample names
        - files, refined samlle's bigwig file paths
    '''
    
    gt = {}
    for k, v in d.items():
        gt[k] = {}
        v = [x for x in v if x in l]
        gt[k]['samples'] = v
        gt[k]['bwfiles'] = [x for x in bwfiles if any([s in x for s in v])]
    
    return gt


def avgBigwig(bwfiles, grange):
    '''Average multiple bigwig files in a specific region
    
    bwfiles : list of bigwig files (path).
    grange  : tuple. Genomic range, BED like 0-based coordinates, 
              eg. ('chr1', 25101, 27101)
    
    return  : a dictionary of keys: 
        - header, for pyBigWig to write as header
        - values, for pyBigwig to addentries as values
    '''
    chrom, start, end = grange
    values = []
    bwo = {}
    for bw in bwfiles:
        if not os.path.isfile(bw): continue
        with pw.open(bw, 'rt') as b:
            header = list(b.chroms().items())
            vals = b.values(chrom, start, end, numpy=True)
            vals = np.nan_to_num(vals)
            values.append(vals)
    
    if values != [] and header != []:
        avgValues = np.mean(values, axis=0)
        bwo = {'header': header, 'values': avgValues}
    return bwo


def writeBigwig(bwo, fout):
    '''Write out bigwig
    bwo : dict. Output of function avgBigwig.
    '''
    



def main(options):
    

    outdir, outbase = os.path.dirname(options.outprefix), os.path.basename(options.outprefix)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    bwfiles = []
    with open(options.bwfiles) as f:
        bwfiles = [x.strip() for x in f.readlines()]

    vcf = ps.VariantFile(options.vcf)
    
    snp = options.snp.strip().split(':')
    snp = snp[0], int(snp[1])
    
    # phenotype range
    a, b = options.region.split(':')
    b, c = b.split('-')
    extend = int(options.extend)
    grange = a, int(b), int(c)
    
    with open(options.samps) as f:
        pSamps = [l.strip() for l in f.readlines() if l.strip() != '']

    # get sample names and bw files
    samples = getSamplesByAllele(vcf, snp)
    samples = refineSamples(samples, pSamps, bwfiles)
    
    # compute average by allele - 0, 1, 2
    avgBW = {}
    for k, d in samples.items():
        chrom, start, end = grange
        avgBW = avgBigwig(d['bwfiles'], (chrom, start-extend, end+extend))
        # bigwig out file
        bwOut = os.path.join(outdir, 
                    '_'.join([outbase, chrom, str(start), str(end), str(k)]) + '.bw'
                            )
        print(bwOut)
        with pw.open(bwOut, 'wt') as bw:
            bw.addHeader(avgBW['header'])
            bw.addEntries(chrom, start - extend, values=avgBW['values'], span=1, step=1)
        
    
    print('\n\n', snp, '\n\n', samples, '\n\n')

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser()
    
    parser.add_argument("-V", "--vcf", dest="vcf",
        type=str, required=True,
        help="gzipped and tbi indexed vcf file with sample genotypes")

    parser.add_argument("-B", "--bwfiles", dest="bwfiles",
        type=str, required=True,
        help="a file that includes a list of Bigwig file paths")
    
    parser.add_argument("-S", "--snp", dest="snp",
        type=str, required=True,
        help="SNP coordinates, following format chr1:1000, 1-based as in vcf.")
    
    parser.add_argument("-P", "--samps", dest="samps",
        type=str, required=True,
        help="A text file with 1 sample ID per line, indicating samples that \
              exist in phenotypes.")
    
    parser.add_argument("-R", "--region", dest="region",
        type=str, required=True,
        help="Phenotype genomic region, format: chr1:1000-2000, 0-base as in bed.")
    
    parser.add_argument("-E", "--extend", dest="extend",
        default = 50, 
        help="Extend region on both ends to help plotting.")
    
    parser.add_argument("-O", "--outprefix", dest="outprefix",
        type=str, required=True,
        help="Output prefix")
    
    options = parser.parse_args()
    
    main(options)