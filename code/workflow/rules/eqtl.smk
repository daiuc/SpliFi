



# ----------------------------------------------------------------------------------------
# run nominal pass eQTL for GTEx
# ----------------------------------------------------------------------------------------

rule PreareGTExEQTLPhenotype:
    input: 
        cnt = 'resources/{datasource}/expression/{group}_gene_reads.tsv.gz',
        genelist = 'resources/{datasource}/ExpressedGeneList.txt'
    output: 
        cpm = 'results/eqtl/{datasource}/{group}/cpm.bed.gz',
        qq = 'results/eqtl/{datasource}/{group}/qqnorm.bed.gz',
    params:
        R_script = 'workflow/scripts/fromBen_PrepareGTExEQTLPheno.R',
    log: 'logs/PrepareGTExEQTLPheno_{datasource}_{group}.log'
    shell:
        '''
        Rscript {params.R_script} {input.cnt} {input.genelist} {output.cpm} {output.qq} 

        '''

rule sortGTExPhenotype:
    input: 'results/eqtl/{datasource}/{group}/qqnorm.bed.gz',
    output: 
        bed = 'results/eqtl/{datasource}/{group}/qqnorm.sorted.bed.gz',
        tbi = 'results/eqtl/{datasource}/{group}/qqnorm.sorted.bed.gz.tbi'
    shell:
        '''
        bedtools sort -header -i {input} | bgzip -c > {output.bed}
        tabix -p bed {output.bed}

        '''
rule GTExPhenoPC:
    input: 'results/eqtl/{datasource}/{group}/qqnorm.sorted.bed.gz'
    output: 'results/eqtl/{datasource}/{group}/qqnorm.sorted.bed.pca'
    params:
        R_script = 'workflow/scripts/fromBen_PermuteAndPCA.R'
    shell:
        '''
        Rscript {params.R_script} {input} {output}
        '''

rule GTEx_eQTL_Nominal:
    message: 'Run nominal pass eQTL for GTEx, with only sQTLs top variants'
    input: 
        vcf = 'results/qtl/noisy/{datasource}/{group}/separateNoise/cis_100000/nom/{chrom}.TopVariants.vcf.gz',
        cov = 'results/eqtl/{datasource}/{group}/qqnorm.sorted.bed.pca',
        bed = 'results/eqtl/{datasource}/{group}/qqnorm.sorted.bed.gz',
    output: temp('results/eqtl/{datasource}/{group}/nom/{chrom}.txt')
    params:
        window = 1000000, # 1Mb
        chrom = '{chrom}'
    log: 'logs/GTEx_eQTL_Nominal_{datasource}_{group}_{chrom}.log'
    shell:
        '''
        module unload gsl && module load gsl/2.5
        QTLtools cis \
            --seed 123 \
            --nominal 1 \
            --vcf {input.vcf} --bed {input.bed} --cov {input.cov}  --out {output} \
            --window {params.window} \
            --region {params.chrom} &> {log}
        '''
use rule TabixNominal as TabixGTExEQtlNominal with:
    input: 'results/eqtl/{datasource}/{group}/nom/{chrom}.txt'
    output: 'results/eqtl/{datasource}/{group}/nom/{chrom}.txt.gz'
    log: 'logs/TabixGTExNominal_{datasource}_{group}_{chrom}.log'

use rule GTEx_eQTL_Nominal as GTEx_eQTL_Nominal_AllSNPS with:
    message: 'Run nominal pass eQTL for GTEx with all genomewide SNPs (cis)'
    input: 
        vcf = 'results/geno/{datasource}/{group}/{chrom}.vcf.gz',
        cov = 'results/eqtl/{datasource}/{group}/qqnorm.sorted.bed.pca',
        bed = 'results/eqtl/{datasource}/{group}/qqnorm.sorted.bed.gz',
    params:
        window = 200000, # 200kb
        chrom = '{chrom}'
    output: temp('results/eqtl/{datasource}/{group}/nom-all-snps/{chrom}.txt')
    log: 'logs/GTEx_eQTL_Nominal_AllSNPS_{datasource}_{group}_{chrom}.log'

use rule TabixNominal as TabixGTExNominal_AllSNPS with:
    input: 'results/eqtl/{datasource}/{group}/nom-all-snps/{chrom}.txt'
    output: 'results/eqtl/{datasource}/{group}/nom-all-snps/{chrom}.txt.gz'
    log: 'logs/TabixGTExNominal_AllSNPS_{datasource}_{group}_{chrom}.log'
