import os
import pandas as pd


configfile: "config/config.yaml"

wildcard_constraints:
    datasource = "Geuvadis|GTEx",
    datasource2 = "Geuvadis|GTEx",
    population = "CEU|FIN|GBR|TSI|YRI|EUR",
    geuvadis_sid = "[A-Z]{2}[0-9]{5,6}",
    chrom = "chr[0-9]{1,2}",
    window = "[0-9]+",
    AB = "\d+_\d+"


regions = pd.read_csv('results/plots/sashimi/Geuvadis/EUR/regions_to_plot.tsv', sep='\t', header=0)
regions.set_index(['range'], drop=False)

CHROMS = regions['chrom']
ABS = regions['range'].str.split('[:-]', expand=False)
ABS = [f'{x[1]}_{x[2]}' for x in ABS]
DATASOURCES = ['Geuvadis'] * len(ABS)
GROUPS = ['EUR'] * len(ABS)

localrules: all

rule all:
    input:
        expand('results/plots/sashimi/{datasource2}/{group}/{chrom}_{AB}/sashimi.pdf', zip,
                datasource2=DATASOURCES, group=GROUPS, chrom=CHROMS, AB=ABS)


rule MakeBigWigFileList:
    message: '### write bigwig file list into file'
    input: 'results/pheno/noisy/{datasource2}/{group}/individuals.txt'
    output: 'results/plots/sashimi/{datasource2}/{group}/bwfiles.txt'
    params:
        prefix = '/project2/yangili1/ankeetashah/hg38_LCL/',
        suffix = '.bam.bw'
    shell:
        '''
        awk '{{print "{params.prefix}"$1"{params.suffix}"}}' {input} > {output}
        '''


def getPlotSashimiVcf(wildcards):
    if wildcards.datasource2 == 'GTEx':
        vcf = config['annotation']['GTEx']['VCF_v7_HG38']
    elif wildcards.datasource2 == 'Geuvadis':
        vcf =  f'/project2/yangili1/zpmu/1kg_b38/CCDG_14151_B01_GRM_WGS_2020-08-05_{wildcards.chrom}.filtered.shapeit2-duohmm-phased.vcf.gz',
    else: 
        print('Error. Invalid vcf file path. Exiting...')
        exit(0)
    return vcf

def getPlotSashimiParams(wildcards):
    '''Shared by both the prep and plot rules
    '''
    ds, grp, chrom, AB = wildcards.datasource2, wildcards.group, wildcards.chrom, wildcards.AB
    A, B = AB.split('_')
    gr = f'{chrom}:{A}-{B}'
    snp = regions.query('range == @gr').reset_index(drop=True)['snp'][0]
    allele = regions.query('range == @gr').reset_index(drop=True)['allele'][0]
    outprefix = f'results/plots/sashimi/{ds}/{grp}/{chrom}_{AB}/plot'

    return {'region': gr, 'snp':snp, 'outprefix':outprefix, 'allele':allele}

rule PrepPlotSashimi:
    message: "Make Sashimi plots"
    input:
        counts = 'results/pheno/noisy/{datasource2}/{group}/leafcutter_perind.counts.noise_by_intron.gz',
        indvs = 'results/pheno/noisy/{datasource2}/{group}/individuals.txt',
        vcf = getPlotSashimiVcf,
        bwlist = 'results/plots/sashimi/{datasource2}/{group}/bwfiles.txt',
        template = 'config/template.ini'
    output:
        outdir = directory('results/plots/sashimi/{datasource2}/{group}/{chrom}_{AB}'),
        done = touch('results/plots/sashimi/{datasource2}/{group}/{chrom}_{AB}/done'),
    params: getPlotSashimiParams
    shell:
        '''
        python workflow/scripts/sashimi_ingredients.py \
            -V {input.vcf} \
            -B {input.bwlist} \
            -C {input.counts} \
            -P {input.indvs} \
            -T {input.template} \
            -S {params[0][snp]} -R {params[0][region]} -E 1000 \
            -O {params[0][outprefix]}
        
        '''

rule PlotSashimi:
    input: 
        flag = 'results/plots/sashimi/{datasource2}/{group}/{chrom}_{AB}/done'
    output: 
        plot = 'results/plots/sashimi/{datasource2}/{group}/{chrom}_{AB}/sashimi.pdf'
    params: 
        d = getPlotSashimiParams,
        ini = 'results/plots/sashimi/{datasource2}/{group}/{chrom}_{AB}/plot.ini'
    conda: 'pygenometracks'
    shell:
        '''
        pyGenomeTracks \
            --tracks {params.ini} \
            --region {params[0][region]} \
            -t "{params[0][region]} ({params[0][allele]})" \
            --width 9 --trackLabelFraction 0.01 \
            -out {output.plot} --fontSize 4

        '''




# rule Sushi:
#     message: 'plot sashimi plots'
#     output: 'mysashimi.pdf'
#     singularity: '/scratch/midway2/chaodai/singularity/ggsashimi.sif'
#     params:
#         bams = 'workflow/submodules/ggsashimi/examples/input_bams.tsv',
#         region = 'chr10:27040584-27048100',
#         gtf = '/project2/yangili1/cdai/genome_index/hs38/gencode.v38.primary_assembly.annotation.gtf.gz',
#         outprefix = 'mysashimi',
#         pyscript = 'workflow/submodules/ggsashimi/ggsashimi.py'
#     shell:
#         '''
#         {params.pyscript} -b {params.bams} -c {params.region} -g {params.gtf} -o {params.outprefix}
#         ls -l {output}
#         '''

# def getPlotSashimiParams(wildcards):
    
#     import pandas as pd
    
#     ds,grp = wildcards.datasource2, wildcards.group
    
#     if wildcards.datasource2 == 'Geuvadis':
#         fin = f'results/plots/sashimi/{ds}/{grp}/regions_to_plot.tsv'
#         df = pd.read_csv(fin, sep='\t', header=0)
#         grange, snp, allele = df['range'], df['snp'], df['allele']
        
#         grange_n = [x.replace(':', '_').replace('-', '_') for x in grange]
#         outprefix = [f'results/plots/sashimi/{ds}/{grp}/{x}/plot' for x in grange_n] 
    
#     return {'grange': grange, 'snp': snp, 'allele': allele, 'outprefix': outprefix}

