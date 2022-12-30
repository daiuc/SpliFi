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
import gzip
import numpy as np
import pysam as ps
import pandas as pd
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

def strTofloat(frac):
    '''Convert string fraction to float
    frac : str. eg. '5/10'
    return : float
    '''
    a, b = frac.split('/')
    a, b = int(a), int(b)
    if b == 0: frac = 0
    else: frac = a/b
    return frac 


def getCountDataframe(countfile, samples, grange, issorted=True):
    '''Get count dataframe for selected samples

    countfile : str. count file, eg. leafcutter_perind.counts.noise_by_intron.gz
    samples   : list. List of selected samples names in phenotype.
    grange    : tuple. eg ('chr1', 1000, 2000)
    sorted    : bool. If True, countfile is sorted. 

    return : a dataframe that contains read counts for selected samples. 
             Following the format required by pygenome track link table.
    '''

    from io import StringIO  
    from statistics import mean
    
    if '.gz' in countfile:
        useGzip = True
        f = gzip.open(countfile)
    else:
        useGzip = False
        f = open(countfile)
    if issorted: print('Count table indicated as sorted.')

    chrom, start, end = grange
        
    i = 0
    stream_N = StringIO() # noisy introns
    stream_F = StringIO() # functional_introns
    for ln in f:
        # if i > 50: break
        if useGzip:
            lnsplit = ln.decode().strip().split()
        else:
            lnsplit = ln.strip().split()
        if i == 0:
            cindx = [getIndex(lnsplit, x)[0] for x in samples] # col index
            oln = ['chr1', 'st1', 'en1', 'chr2', 'st2', 'en2', 'meanPSI'] # line out
            buf = '\t'.join(oln) + '\n'
            stream_N.write(buf)
            stream_F.write(buf)
        else:
            if chrom == lnsplit[0].split(":")[0]: # match chromosome
                oln = [lnsplit[j] for j in [0] + cindx] # line out
                ch, st, en, cl, an = lnsplit[0].split(':')
                _, _, stra = cl.split('_')

                if int(st) >= start and int(en) <= end: # match range
                    meanPSI = mean([strTofloat(x) for x in oln[1:]])
                    # oln = [ch, st, en, lnsplit[0], stra, an, str(meanPSI)]
                    oln = [ch, st, st, ch, en, en, str(meanPSI)]
                    buf = '\t'.join(oln) + '\n'
                    if an == "N": stream_N.write(buf)
                    elif an == "F": stream_F.write(buf)
                elif issorted and int(st) > end+10: # line passed range end loop
                    break
                else: continue
            else: continue
        i += 1
    
    stream_N.seek(0)
    stream_F.seek(0)
    df_N = pd.read_csv(stream_N, sep='\t')
    df_F = pd.read_csv(stream_F, sep='\t')
    stream_N.close()
    stream_F.close()
    return {'N':df_N[df_N.meanPSI > 0, 'F':df_F[df_F.meanPSI > 0]


def writeIni(template, outfile, fdict):
    '''Generate ini file for plotting sashimi plots
    template : str. Path to a template ini file.
    outfile  : str. Path of ini file to write into.
    fdict.   : dict. Dictionary hosting bw file path, and link file path
                eg {0: {'bw': 'allele0.bw', 'link': 'lnk0.tsv'}, 1:...}
    
    return : None.
    '''
    
    from jinja2 import Template
    
    print(f'Generate ini files: {outfile} ...\n')
    with open(template) as t:
        t_content = t.read()
        tt = Template(t_content)
        outini = tt.render( # out ini
                title = "sashimi",
                title_height = 1,
                where = "right",
                title_fontsize = 4,
                allele0 = os.path.basename(fdict[0]['bw']).split('.')[1],
                bw0 = fdict[0]['bw'],
                lnk0 = fdict[0]['link'],
                allele1 = os.path.basename(fdict[1]['bw']).split('.')[1],
                bw1 = fdict[1]['bw'],
                lnk1 = fdict[1]['link'],
                allele2 = os.path.basename(fdict[2]['bw']).split('.')[1],
                bw2 = fdict[2]['bw'],
                lnk2 = fdict[2]['link'],
        )
    with open(outfile, 'w') as f:
        f.write(outini + '\n')

def main(options):
    
    outdir, outbase = os.path.dirname(options.outprefix), os.path.basename(options.outprefix)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    # bigwig files
    bwfiles = []
    with open(options.bwfiles) as f:
        bwfiles = [x.strip() for x in f.readlines()]

    # indexed bcf file for extracting genotypes
    vcf = ps.VariantFile(options.vcf)
    
    # SNP
    snp = options.snp.strip().split(':')
    snp = snp[0], int(snp[1])
    
    # phenotype range, use extend to extend both ends
    a, b = options.region.split(':')
    b, c = b.split('-')
    extend = int(options.extend)
    grange = a, int(b), int(c)
    
    # sample names in both phenotypes and genotypes
    with open(options.samps) as f:
        pSamps = [l.strip() for l in f.readlines() if l.strip() != '']

    # get sample names and bw files
    samples = getSamplesByAllele(vcf, snp)
    samples = refineSamples(samples, pSamps, bwfiles)
    
    
    avgBW = {} # avg bigwig for an allele
    outDict = {} # output file for each allele 0, 1, 2
    for k, d in samples.items():
        chrom, start, end = grange
        
        # Average bigwig by allele and write average bw files
        avgBW = avgBigwig(d['bwfiles'], (chrom, start-extend, end+extend))
        bwOut = os.path.join(outdir, f'{outbase}_{str(k)}.bw')

        print(f'Extracting average bigwig for allele {k}, write to {bwOut}...\n')
        with pw.open(bwOut, 'wt') as bw:
            bw.addHeader(avgBW['header'])
            bw.addEntries(chrom, start - extend, values=avgBW['values'], span=1, step=1)
        
        # Extract PSI for allele and write link table
        lnkOut_N = os.path.join(outdir, f'{outbase}_Noisy_{str(k)}.tsv')
        lnkOut_F = os.path.join(outdir, f'{outbase}_Func_{str(k)}.tsv')

        print(f'Extract link tables for allele {k}, write to {lnkOut_N}, {lnkOut_F}...\n')

        dfs = getCountDataframe(options.countfile, d['samples'], 
                (chrom, start-extend, end+extend), options.issorted)
        
        if not dfs['N'].empty:
            dfs['N'].to_csv(lnkOut_N, sep="\t", header=False, index=False)
        else:
            print(f"Allele {k} doesn't have noisy introns with any counts within selected region.\n")

        if not dfs['F'].empty:
            dfs['F'].to_csv(lnkOut_F, sep="\t", header=False, index=False)
        else:
            print(f"Allele {k} doesn't have functional introns with any counts within selected region.\n")

        # store written bw and link file paths by allele for ini generation, note only basename for bw file in ini
        outDict[k] = {'bw': os.path.basename(bwOut), 
                      'link_N': os.path.basename(lnkOut_N),
                      'link_F': os.path.basename(lnkOut_F),
                      }
        
        iniOut = os.path.join(outdir, f'{outbase}.ini')

    # Generate ini file for sashimi plot
    if options.template:
        writeIni(options.template, iniOut, outDict)

    # print('\n\n', snp, '\n\n', samples, '\n\n')
    # print(linkdf)

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser()
    
    parser.add_argument("-V", "--vcf", dest="vcf",
        type=str, required=True,
        help="gzipped and tbi indexed vcf file with sample genotypes")

    parser.add_argument("-B", "--bwfiles", dest="bwfiles",
        type=str, required=True,
        help="a file that includes a list of Bigwig file paths")
    
    parser.add_argument("-C", "--countfile", dest="countfile",
        type=str, required=True, help="a count matrix file")
    
    parser.add_argument("-D", "--issorted", dest="issorted", default=True,
        action="store_false", help="Sorted count file speeds things up. (default: true)")


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
    
    parser.add_argument("-T", "--template", dest="template",
        type=str, required=False, help="ini template")

    parser.add_argument("-O", "--outprefix", dest="outprefix",
        type=str, required=True, help="Output prefix")
    
    options = parser.parse_args()
    
    main(options)