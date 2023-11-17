'''
Rules for running QTL analysis

NOTE qtltools require gsl, make sure gsl is loaded before running snakemake
'''

#------------------------------------------------------------------#
#             QTLtools using *.counts.noise.gz phenotype           #
#------------------------------------------------------------------#


def getExtractNoisyPhenoInput(wildcards):

    counts = os.path.join('results/pheno/noisy', wildcards.datasource, 
                          wildcards.group, 'leafcutter_perind.counts.noise.gz')
    counts_by_intron = os.path.join('results/pheno/noisy', wildcards.datasource, 
                          wildcards.group, 'leafcutter_perind.counts.noise_by_intron.gz')
    if wildcards.datasource == 'GTEx':
        vcf = config['annotation']['GTEx']['VCF_v7_HG38']
    elif wildcards.datasource == 'Geuvadis':
        # only need the vcf header to get all available sample names, thus choose a smaller vcf
        vcf = '/project/yangili1/zpmu/1kg_b38/CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz'
    return {'counts': counts, 'counts_by_intron': counts_by_intron, 'vcf': vcf}


rule ExtractNoisyCounts:
    input: unpack(getExtractNoisyPhenoInput)
    output:
        counts_noise = 'results/pheno/noisy/{datasource}/{group}/leafcutter_perind.extract.agg_noise_only.txt.gz',
        counts_both = 'results/pheno/noisy/{datasource}/{group}/leafcutter_perind.extract.agg_both.txt.gz',
        counts_by_intron = 'results/pheno/noisy/{datasource}/{group}/leafcutter_perind.extract.by_intron.txt.gz',
        indivs = 'results/pheno/noisy/{datasource}/{group}/individuals.txt'
    params:
        py_script = 'workflow/scripts/extract_leafcutter_counts_by_vcf.py',
        out_prefix = 'results/pheno/noisy/{datasource}/{group}/leafcutter_perind'
    shell:
        '''
        python {params.py_script} --VCF {input.vcf} --counts {input.counts} \
            --counts_by_intron {input.counts_by_intron} \
            --OutPrefix {params.out_prefix} --OutIndivs {output.indivs}

        ls {output.counts_noise} {output.counts_both} {output.counts_by_intron} {output.indivs}
    
        '''


rule PrepPhenoBed:
    message: '### Prepare phenotype bed file for qtltools (rank normalized)'
    input:
        counts = rules.ExtractNoisyCounts.output.counts_both,
        anno = config['annotation']['gencode']
    # output: 'results/pheno/noisy/{datasource}/{group}/pheno.chr1.bed.gz'
    output: touch('results/pheno/noisy/{datasource}/{group}/PrepPhenoBed.done')
    params:
        Rscript = 'workflow/scripts/prepPhenoBed.R',
        outprefix = 'results/pheno/noisy/{datasource}/{group}/pheno',
        minclu = 10, # min mean cluster reads
        minsam = 10, # min number of samples passing
        minread = 2  # min number of reads per sample
    threads: 8
    resources: cpu = 8, mem_mb = 28000, time = '12:00:00'
    shell:
        '''
        Rscript {params.Rscript} -I {input.counts} -A {input.anno} -O {params.outprefix} \
            -C {params.minclu} -S {params.minsam} -N {params.minread} -T {threads}
        
        beds=({params.outprefix}.*.bed)
        for b in ${{beds[@]}}; do
            bgzip -f $b
            tabix -p bed ${{b}}.gz
        done
        
        '''



rule PhenotypePCA:
    message: '### Run PCA on phenotype with permutation'
    input: rules.PrepPhenoBed.output
    output: 'results/pheno/noisy/{datasource}/{group}/{chrom}.pca'
    params:
        rscript = 'workflow/scripts/PermuteAndPCA.R',
        inputfile = 'results/pheno/noisy/{datasource}/{group}/pheno.{chrom}.bed.gz'
    shell:
        '''
        ls {input}
        Rscript {params.rscript} {params.inputfile} {output} 
        '''



def getExtractGenotypeInput(wildcards):
    if wildcards.datasource == 'GTEx':
        vcf = config['annotation']['GTEx']['VCF_v7_HG38']
    elif wildcards.datasource == 'Geuvadis':
        vcf =  f'/project/yangili1/zpmu/1kg_b38/CCDG_14151_B01_GRM_WGS_2020-08-05_{wildcards.chrom}.filtered.shapeit2-duohmm-phased.vcf.gz',
    else: 
        print('Error. Invalid vcf file path. Exiting...')
        exit(0)
    return vcf

def getExtractGenotypeParams(wildcards):
    if wildcards.datasource == 'GTEx':
        min_MAF, chrom = 0.05, wildcards.chrom
        expr = f'AF >= {min_MAF} ' 
    elif wildcards.datasource == 'Geuvadis':
        min_MAF, max_HWE = 0.05, 1e-3
        expr = f'AF >= {min_MAF} && HWE < {max_HWE}' 
    else: 
        print('Error. Invalid vcf file path. Exiting...')
        exit(0)
    return expr

rule ExtractGenotypeVCF:
    '''
        Procedures:
            1.  extract vcf per chromosome for samples in phenotype
            2.  tabix vcf
    '''
    message: '### Extract genotype vcf for phenotype samples - {wildcards.datasource}:{wildcards.group}:{wildcards.chrom}'
    input: 
        vcf = getExtractGenotypeInput,
        sample_file = rules.ExtractNoisyCounts.output.indivs
    output: 'results/geno/{datasource}/{group}/{chrom}.vcf.gz'
    params:
        expr = getExtractGenotypeParams
    threads: 4
    resources: time=2000, mem_mb=15000, cpu=4
    shell:
        '''
        bcftools view \
            --threads {threads} \
            --samples-file {input.sample_file} \
            -i "{params.expr}" \
            -Oz -o {output} \
            {input.vcf} \
            {wildcards.chrom}
        
        bcftools index --threads {threads} --tbi {output}
        '''

# sometimes qtltools fail, module unload gsl gcc and reload gcc/10.2.0 and gsl/2.2.1 works
rule GenotypePCA:
    input:  'results/geno/{datasource}/{group}/{chrom}.vcf.gz'
    output: 'results/geno/{datasource}/{group}/{chrom}.pca'
    params:
        out_prefix = 'results/geno/{datasource}/{group}/{chrom}'
    shell:
        '''
            QTLtools pca \
            --seed 123 \
            --maf 0.01 \
            --vcf {input}  \
            --out {params.out_prefix} \
            --center \
            --scale
        '''


rule MakeCovarianceMatrix:
    message: '### Make covariance matrix with fixed PCs'
    input:
        PhenoPCs = rules.PhenotypePCA.output,
        GenoPCs = rules.GenotypePCA.output
    output: 'results/pheno/noisy/{datasource}/{group}/{chrom}_CovMatrix.txt'
    params:
        Geno_PCs = 4,
        rscript = 'workflow/scripts/Make_CovMatrix.R'
    shell:
        '''
        Rscript {params.rscript} {input.PhenoPCs} {input.GenoPCs} {params.Geno_PCs} {output}
        '''


rule MapQTL_Perm:
    '''
        already ranknorm normalized, do not need the normal option 
    '''
    message: 'Map QTL using permutation pass'
    input:
        flag = 'results/pheno/noisy/{datasource}/{group}/PrepPhenoBed.done',
        vcf = 'results/geno/{datasource}/{group}/{chrom}.vcf.gz',
        cov = 'results/pheno/noisy/{datasource}/{group}/{chrom}_CovMatrix.txt'
    output: temp('results/qtl/noisy/{datasource}/{group}/cis_{window}/perm/{chrom}.txt')
    params:
        bed = 'results/pheno/noisy/{datasource}/{group}/pheno.{chrom}.bed.gz',
        cis_window = '{window}'
    resources: cpu = 1, mem = 12000, time = 1000
    shell:
        '''
        ls {input.flag}
        QTLtools cis \
            --seed 123 --silent \
            --vcf {input.vcf} --bed {params.bed} --cov {input.cov}  --out {output} \
            --window {params.cis_window} \
            --permute 1000 --region {wildcards.chrom}
        '''


rule AddQvalueToPermutationPass:
    message: '### Add qvalue to the output of QTLtools permutation pass'
    input:  'results/qtl/noisy/{datasource}/{group}/cis_{window}/perm/{chrom}.txt'
    output: 'results/qtl/noisy/{datasource}/{group}/cis_{window}/perm/{chrom}.addQval.txt.gz'
    params:
        rscript = 'workflow/scripts/AddQvalueToQTLtoolsOutput.R'
    shell:
        '''
            txt=$(echo {output} | sed -E 's/.gz//')
            Rscript {params.rscript} {input} $txt
            bgzip $txt
        '''

    # QTLtools cis permutation pass output fields:
    # - 1. phenotype_id
    # - 2. phenotype_chr
    # - 3. phenotype_start
    # - 4. phenotype_end
    # - 5. phenotype_strand
    # - 6. number_variants_tested
    # - 7. best_nominal_distance
    # - 8. best_genotype_id
    # - 9. best_genotype_chr
    # - 10. best_genotype_start
    # - 11. best_genotype_end
    # - 12. true_degrees_of_freedom
    # - 13. estimated_dof
    # - 14. beta_ml1
    # - 15. beta_ml2
    # - 16. pval_nominal
    # - 17. pval_r2
    # - 18. regression_slop
    # - 19. pval_empirical (by direct permutation method)
    # - 20. pval_adjusted (adjusted p-value, note, this has not been genome-wide multiple-tested on phenotypes)




#------------------------------------------------------------------#
#        Map QTLs using *.counts.noise_by_intron.gz phenotype      #
#------------------------------------------------------------------#


use rule PrepPhenoBed as prepPhenoBed_by_intron with:
    input:
        counts = rules.ExtractNoisyCounts.output.counts_by_intron,
        anno = config['annotation']['gencode']
    output: touch('results/pheno/noisy/{datasource}/{group}/by_intron/PrepPhenoBed.done')
    params:
        Rscript = 'workflow/scripts/prepPhenoBed.R',
        outprefix = 'results/pheno/noisy/{datasource}/{group}/by_intron/pheno',
        minclu = 10, # min mean cluster reads
        minsam = 10, # min number of samples passing
        minread = 2  # min number of reads per sample
    threads: 8
    resources: cpu = 8, mem_mb = 28000, time = '12:00:00'


use rule PhenotypePCA as PhenotypePCA_by_intron with:
    input: rules.prepPhenoBed_by_intron.output
    output: 'results/pheno/noisy/{datasource}/{group}/by_intron/{chrom}.pca'
    params:
        rscript = 'workflow/scripts/PermuteAndPCA.R',
        inputfile = 'results/pheno/noisy/{datasource}/{group}/by_intron/pheno.{chrom}.bed.gz'


use rule MakeCovarianceMatrix as MakeCovarianceMatrix_by_intron with:
    input:
        PhenoPCs = rules.PhenotypePCA_by_intron.output,
        GenoPCs = rules.GenotypePCA.output # use the same vcf as before
    output: 'results/pheno/noisy/{datasource}/{group}/by_intron/{chrom}_CovMatrix.txt'
    params:
        Geno_PCs = 4,
        rscript = 'workflow/scripts/Make_CovMatrix.R'


use rule MapQTL_Perm as MapQTL_Perm_by_intron with:
    input:
        flag = 'results/pheno/noisy/{datasource}/{group}/by_intron/PrepPhenoBed.done',
        vcf = 'results/geno/{datasource}/{group}/{chrom}.vcf.gz',
        cov = 'results/pheno/noisy/{datasource}/{group}/by_intron/{chrom}_CovMatrix.txt'
    output: temp('results/qtl/noisy/{datasource}/{group}/by_intron/cis_{window}/perm/{chrom}.txt')
    params:
        bed = 'results/pheno/noisy/{datasource}/{group}/by_intron/pheno.{chrom}.bed.gz',
        cis_window = '{window}'
    resources: cpu = 1, mem = 12000, time = 1000


use rule AddQvalueToPermutationPass as AddQvalueToPermutationPass_by_intron with:
    input: rules.MapQTL_Perm_by_intron.output
    output: 'results/qtl/noisy/{datasource}/{group}/by_intron/cis_{window}/perm/{chrom}.addQval.txt.gz'
    params:
        rscript = 'workflow/scripts/AddQvalueToQTLtoolsOutput.R'



