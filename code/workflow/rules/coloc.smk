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

# this rule is used to subset expression and splicing phenotype ids for nominal pass run
# It select phenotypes that have a significant eQTL/sQTL; and for splicing, the phenotype
# must be either PR or UP or NE.
# It output a table listing selected splicing pid and expression pid, as well as
# qqnormed expression and splicing phenotype files, ready for nominal pass run.
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


#--- Nominal pass with SE, output columns ---
# col1 pid
# col2 pchr
# col3 pstart
# col4 pend
# col5 pstrand
# col6 number of snps
# col7 distance
# col8 snp id
# col9 snp chr
# col10 snp start
# col11 snp end
# col12 nominal pval
# col13 r2
# col14 slope
# col15 slope SE
# col16 1/0 to indicate top snp


# this rule run nominal pass for expression on selected genes
# The output provides summary stats for coloc analysis.
rule eQTLNominalForColoc:
    message: 'Map eQTL nominal pass for eQTL-sQTL coloc, report beta SE'
    input: 
        vcf = 'results/geno/{datasource}/{group}/{chrom}.vcf.gz',
        cov = 'results/eqtl/{datasource}/{group}/qqnorm.sorted.{chrom}.cov',
        bed = 'results/coloc/sqtl-eqtl/{datasource}/{group}/{chrom}.eqtl_pheno.bed.gz',
    output: 'results/coloc/sqtl-eqtl/{datasource}/{group}/{chrom}.eqtl_nominal.txt.gz'
    params:
        window = 1500000, # 1.5Mb
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
            --window {params.window} --silent

        awk 'BEGIN {{OFS="\t"}} {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}}' $(echo {output} | sed -e 's/\.gz//') | \
            sort -k9,9V -k10,10n -k11,11n | \
            bgzip -c > {output}

        tabix -s 9 -b 10 -e 11 -f {output}
        '''

# same as the above rule, but for splicing
use rule eQTLNominalForColoc as sQTLNominalForColoc with:
    message: 'Map sQTL nominal pass for eQTL-sQTL coloc, report beta SE'
    input: 
        vcf = 'results/geno/{datasource}/{group}/{chrom}.vcf.gz',
        cov = 'results/pheno/noisy/{datasource}/{group}/separateNoise/{chrom}_CovMatrix.txt',
        bed = 'results/coloc/sqtl-eqtl/{datasource}/{group}/{chrom}.sqtl_pheno.bed.gz',
    output: 'results/coloc/sqtl-eqtl/{datasource}/{group}/{chrom}.sqtl_nominal.txt.gz'
    params:
        window = 1500000, # also use 1.5Mb window
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
        coloc_res_full = 'results/coloc/sqtl-eqtl/{datasource}/{group}/{chrom}.coloc-res.rds',
    resources: 
        cpu = 1, mem_mb = 25000, time = 1500
    shell:
        '''
        Rscript {params.Rscript} \
            {input.coloc_eqtl_sqtl_ids} \
            {input.eqtl_nominal} \
            {input.sqtl_nominal} \
            {params.coloc_res} \
            {params.coloc_res_full}
        '''


# GWAS hit loic and munged gwas stats are processed using '.../splice-pub/analysis/2024-08-06-process-GWAS-summary-stats.qmd'

rule MakeGwasHitLoci:
    input: 
        gwas_stats = 'resources/GWAS/hg38/{gtrait}.tsv.gz'
    output: 'resources/GWAS/hg38/{gtrait}_hitloci_1en5.bed' # 1en5 = 1e-5
    params:
        Rscript = 'workflow/scripts/GwasHits-v3.R',
        MinPval = 1e-5,
    resources:
        cpu = 1, mem_mb = 25000, time = 1500
    log: 'logs/gwas/{gtrait}.log'
    shell:
        '''
        Rscript {params.Rscript} {input} {output} {params.MinPval} &> {log}
        '''


rule HyPrColoc_Gwas_eQTL_sQTL:
    input:
        gwas_loci = 'resources/GWAS/hg38/{gtrait}_hitloci_1en5.bed',
        gwas_stats = 'resources/GWAS/hg38/{gtrait}.tsv.gz'
    output: touch('results/coloc/gwas-eqtl-sqtl/{datasource}/{gtrait}/{group}/done')
    params:
        gene_annot = '/project/yangili1/cdai/annotations/hg38/gencode.v26.GRCh38.genes.csv',
        eqtl_sqtl_prefix = 'results/coloc/sqtl-eqtl/{datasource}/{group}',
        rscript = 'workflow/scripts/HyprColoc_gwas_eQTL_sQTL.R',
        out_dir = 'results/coloc/gwas-eqtl-sqtl/{datasource}/{gtrait}/{group}',
    log: 'logs/coloc/gwas-eqtl-sqtl/{datasource}/{gtrait}/{group}.log'
    resources: 
        cpu = 1, mem_mb = 25000, time = 2100
    shell:
        '''
        Rscript {params.rscript} \
            {input.gwas_loci} \
            {input.gwas_stats} \
            {params.gene_annot} \
            {params.eqtl_sqtl_prefix} \
            {params.out_dir} &> {log}
        '''


rule computeLD:
    input:
        vcf = config['VCF']['GTEx']['HG38_v7']
        # vcf = "results/geno/GTEx/Lung/chr22.vcf.gz"
    output: 
        bed = "resources/GTEx/LD/gtex.bed"
    params:
        bed_prefix = "resources/GTEx/LD/gtex",
        ld_prefix = "resources/GTEx/LD/gtex_ld_matrix",
    threads: 8
    resources:
        cpu = 8, mem_mb = 35000, time = 2100
    log: 'logs/LD/gtex.log'
    shell:
        '''
        module load plink

        plink --make-bed --vcf {input.vcf} --out {params.bed_prefix}
        sleep 5
        plink --bfile {params.bed_prefix} --r2 --ld-window-r2 0.2 --out {params.ld_prefix} --threads {threads}


        '''






















