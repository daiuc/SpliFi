#!/bin/env python

import argparse
import gzip
import subprocess


def args_parser():
    
    
    parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--VCF",
            help= "gzipped and tbi indexed vcf file with sample genotypes",
            metavar="<VCF file name>",
            default=None)
    parser.add_argument("--counts",
            help= "leafcutter output *.counts.noise.gz file",
            metavar="<counts-file>",
            default=None)
    parser.add_argument("--counts_by_intron",
            help= "leafcutter output *.counts.noise_by_intron.gz file",
            metavar="<counts_by_intron-file>",
            default=None)
    parser.add_argument("--OutPrefix",
            help= "Output prefix for all output files",
            metavar="OutputPrefix",
            default=None)
    parser.add_argument("--OutIndivs",
            help= "Output file with individual ids",
            metavar="OutputIndivs",
            default=None)

    return parser.parse_args()


def main(args):
    
    vcf_file = args.VCF
    counts_file = args.counts
    counts_intron_file = args.counts_by_intron
    out_prefix = args.OutPrefix

    # extract individual ids from vcf file
    vcf_header = f"bcftools view -h {vcf_file} | awk '$1 ~ /#CHROM/' "
    vcf_header = subprocess.run(vcf_header, shell=True, capture_output=True).stdout.decode('utf-8').strip().split()
    geno_inds = vcf_header[9:] # available genotype individual IDs

    # output files
    agg_noise_only = f'{out_prefix}.extract.agg_noise_only.txt.gz'
    agg_both = f'{out_prefix}.extract.agg_both.txt.gz'
    by_intron = f'{out_prefix}.extract.by_intron.txt.gz'
    indivs = args.OutIndivs


    fout_agg_noise_only = gzip.open(agg_noise_only, 'wt')
    fout_agg_both = gzip.open(agg_both, 'wt')
    fout_by_intron = gzip.open(by_intron, 'wt')
    
    # extract counts for aggregated productive/unproductive introns
    i = 0
    clus = []
    for ln in gzip.open(counts_file): # counts.noise.gz
        lnsplit = ln.decode().strip().split()
        if i == 0: # header
            keep1 = [lnsplit.index(l) for l in lnsplit if l in ['chrom'] + geno_inds ] # cols to keep
            keep2 = [lnsplit.index(l) for l in lnsplit if l in geno_inds ] # cols to keep, no 'chrom'
            if len(keep2) < 3:
                print(f'Error. Only {len(keep2)} samples left after removing samples not in genotype vcf. Exiting...')
                exit(0)
            with open(indivs, 'w') as sf: # write individual ids
                indvs = [lnsplit[ix] for ix in keep2]
                print(f'Removing these samples not in genotype vcf:\n{[x for x in lnsplit[1:] if x not in indvs]} ...')
                sf.writelines([x + '\n' for x in indvs])
            buf = ' '.join([lnsplit[ix] for ix in keep1]) + '\n' # header line
            fout_agg_noise_only.write(buf)
            fout_agg_both.write(buf)
        elif '*' in lnsplit[0]:
            # output unproductive counts (aggregated by cluster)
            buf = ' '.join([lnsplit[ix] for ix in keep1]) + '\n'
            # only select clusters that have both productive and unproductive introns
            # extract cluster id, eg. 'chr1:939412:942409:clu_6_+_*' to 'clu_6_+'
            clu = lnsplit[0].split(':')[3].split('_')[0:3]
            clu = '_'.join(clu)
            clus.append(clu)
            fout_agg_noise_only.write(buf)
            fout_agg_both.write(buf)
        else:
            buf = ' '.join([lnsplit[ix] for ix in keep1]) + '\n'
            fout_agg_both.write(buf)
        i += 1
        
    # extract counts for individual productive/unproductive introns
    i = 0
    for ln in gzip.open(counts_intron_file): # counts.noise_by_intron.gz
        lnsplit = ln.decode().strip().split()
        if i == 0: # header
            keep1 = [lnsplit.index(l) for l in lnsplit if l in ['chrom'] + geno_inds ] # cols to keep
            buf = ' '.join([lnsplit[ix] for ix in keep1]) + '\n' # header line
            fout_by_intron.write(buf)
        else:
            clu = lnsplit[0].split(':')[3]
            # if clu in clus: # only clusters that have both productive and unproductive introns
            #     buf = ' '.join([lnsplit[ix] for ix in keep1]) + '\n'
            #     fout_by_intron.write(buf)
            buf = ' '.join([lnsplit[ix] for ix in keep1]) + '\n'
            fout_by_intron.write(buf)
        i += 1

    fout_agg_noise_only.close()
    fout_agg_both.close()
    fout_by_intron.close()


if __name__ == "__main__":
    
    args = args_parser()
    main(args)



