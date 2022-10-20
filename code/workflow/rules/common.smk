import pandas as pd
import os


ERRs = []
with open('resources/geuvadis-subset79-ERRs.txt') as f:
    ERRs = f.readlines()
ERRs = [x.strip() for x in ERRs]

CHROMS = ['chr'+str(i) for i in range(1,23)]


rule Get1KGP_SampleIDs:
    message: 'Collect all Sample IDs that exist in 1000 Genome Project Phase 3 VCF file'
    input: 
        vcf = '/project2/yangili1/zpmu/1kg_b38/CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.filtered.shapeit2-duohmm-phased.vcf.gz'
    output: 
        'resources/SampleID_List_1000GenomePhase3.txt'
    threads: 2
    resources: cpu = 2, mem_mb = 15000, time = 2100
    shell: 
        '''
        bcftools view -h {input.vcf} | \
            awk ''
        '''