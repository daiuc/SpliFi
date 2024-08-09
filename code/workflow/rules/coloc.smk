'''
coloc with eQTLs, sQTLs, GWAS

'''


rule MungeGwasSummaryStats:
    '''Interactive ONLY!'''
    output: 'resources/GWAS/MungedSumStats/done'
    params:
        gwas_yaml = 'resources/GWAS/selected-gwas.yaml',
        min_p_val = 1e-5,
        out_dir = 'resources/GWAS/MungedSumStats',
    notebook: '../scripts/MungeGwasSummaryStats.r.ipynb'


rule SqtlEqtlPhenoForColoc:
    message: 'Prepare sQTL and eQTL data for nominal pass for coloc'
    input:
        eqtl_perm = 'results/eqtl/{datasource}/{group}/perm/{chrom}.txt.addQval.txt.gz',
        sqtl_perm = 'results/qtl/noisy/{datasource}/{group}/separateNoise/cis_100000/perm/{chrom}.addQval.txt.gz',
        eqtl_pheno = 'results/eqtl/{datasource}/{group}/qqnorm.sorted.{chrom}.bed.gz',
        sqtl_pheno = 'results/pheno/noisy/{datasource}/{group}/separateNoise/leafcutter.qqnorm_{chrom}.gz',
    output:
        coloc_eqtl_sqtl_ids = 'results/coloc/sqtl-eqtl/{datasource}/{group}/{chrom}.eqtl_sqtl_ids.txt',
        eqtl_pheno = 'results/coloc/sqtl-eqtl/{datasource}/{group}/{chrom}.eqtl_pheno.bed.gz',
        sqtl_pheno = 'results/coloc/sqtl-eqtl/{datasource}/{group}/{chrom}.sqtl_pheno.bed.gz',
    params:
        gene_annot = '/project/yangili1/cdai/annotations/hg38/gencode.v26.GRCh38.genes.csv',
        fdr = 0.1,
        Rscript = 'workflow/scripts/prepareColoc.R',
    shell:
        '''
        Rscript {params.Rscript} \
            {input.eqtl_perm} \
            {input.sqtl_perm} \
            {params.fdr} \
            {params.gene_annot} \
            {input.eqtl_pheno} \
            {input.sqtl_pheno} \
            $(echo {output.eqtl_pheno} | sed -e 's/\.gz//') \
            $(echo {output.sqtl_pheno} | sed -e 's/\.gz//') \
            {output.coloc_eqtl_sqtl_ids}

        
        bgzip $(echo {output.eqtl_pheno} | sed -e 's/\.gz//')
        bgzip $(echo {output.sqtl_pheno} | sed -e 's/\.gz//')

        tabix -p bed {output.eqtl_pheno}
        tabix -p bed {output.sqtl_pheno}

        '''

rule eQTLNominalForColoc:
    message: 'Map eQTL nominal pass for eQTL-sQTL coloc, report beta SE'
    input: 
        vcf = 'results/geno/{datasource}/{group}/{chrom}.vcf.gz',
        cov = 'results/eqtl/{datasource}/{group}/qqnorm.sorted.{chrom}.cov',
        bed = 'results/coloc/sqtl-eqtl/{datasource}/{group}/{chrom}.eqtl_pheno.bed.gz',
    output: 'results/coloc/sqtl-eqtl/{datasource}/{group}/{chrom}.eqtl_nominal.txt.gz'
    params:
        window = 1000000, # 1Mb
        Rscript = 'workflow/scripts/FormatQtltoolsNominalResult.R',
    log: 'logs/coloc/sqtl-eqtl/{datasource}/{group}/{chrom}.eqtl_nominal.log'
    shell:
        '''
        module unload gsl && module load gsl/2.5

        QTLtools cis \
            --seed 123 \
            --nominal 1 --std-err \
            --vcf {input.vcf} --bed {input.bed} --cov {input.cov}  \
            --out $(echo {output} | sed -e 's/\.gz//') \
            --window {params.window}

        bgzip $(echo {output} | sed -e 's/\.gz//')
        tabix -p bed {output}
        '''

use rule eQTLNominalForColoc as sQTLNominalForColoc with:
    message: 'Map sQTL nominal pass for eQTL-sQTL coloc, report beta SE'
    input: 
        vcf = 'results/geno/{datasource}/{group}/{chrom}.vcf.gz',
        cov = 'results/pheno/noisy/{datasource}/{group}/separateNoise/{chrom}_CovMatrix.txt',
        bed = 'results/coloc/sqtl-eqtl/{datasource}/{group}/{chrom}.sqtl_pheno.bed.gz',
    output: 'results/coloc/sqtl-eqtl/{datasource}/{group}/{chrom}.sqtl_nominal.txt.gz'
    params:
        window = 1000000, # also use 1Mb window
    log: 'logs/coloc/sqtl-eqtl/{datasource}/{group}/{chrom}.sqtl_nominal.log'


rule HyPrColoc_eQTL_sQTL:
    input:
        eqtl_nominal = 'results/coloc/sqtl-eqtl/{datasource}/{group}/{chrom}.eqtl_nominal.txt.gz',
        sqtl_nominal = 'results/coloc/sqtl-eqtl/{datasource}/{group}/{chrom}.sqtl_nominal.txt.gz',
        coloc_eqtl_sqtl_ids = 'results/coloc/sqtl-eqtl/{datasource}/{group}/{chrom}.eqtl_sqtl_ids.txt',
    output: touch('results/coloc/sqtl-eqtl/{datasource}/{group}/{chrom}.coloc.done')
    params:
        Rscript = 'workflow/scripts/HyprColoc_eQTL_sQTL.R',
        coloc_res = 'results/coloc/sqtl-eqtl/{datasource}/{group}/{chrom}.coloc-res.txt',
    log: 'logs/coloc/sqtl-eqtl/{datasource}/{group}/{chrom}.coloc.log'
    shell:
        '''
        Rscript {params.Rscript} \
            {input.coloc_eqtl_sqtl_ids} \
            {input.eqtl_nominal} \
            {input.sqtl_nominal} \
            {params.coloc_res}
        '''























