



# ----------------------------------------------------------------------------------------
# qtltools input prep
# ----------------------------------------------------------------------------------------

rule PreareGTExEQTLPhenotype:
    input: 
        cnt = 'resources/{datasource}/expression/{group}_gene_reads.tsv.gz',
        genelist = 'resources/{datasource}/ExpressedGeneList.txt',
        samples = 'results/pheno/noisy/{datasource}/{group}/separateNoise/leafcutter_names.txt' # intersect of gtex pheno and gtex vcf avail. samples
    output: 
        cpm = 'results/eqtl/{datasource}/{group}/cpm.bed.gz',
        qq = 'results/eqtl/{datasource}/{group}/qqnorm.bed.gz',
    params:
        R_script = 'workflow/scripts/fromBen_PrepareGTExEQTLPheno.R',
    log: 'logs/PrepareGTExEQTLPheno_{datasource}_{group}.log'
    shell:
        '''
        Rscript {params.R_script} {input.cnt} {input.genelist} {input.samples} {output.cpm} {output.qq} 

        '''

rule sortGTExPhenotype:
    input: 'results/eqtl/{datasource}/{group}/qqnorm.bed.gz',
    output: 
        bed = 'results/eqtl/{datasource}/{group}/qqnorm.sorted.bed.gz',
        tbi = 'results/eqtl/{datasource}/{group}/qqnorm.sorted.bed.gz.tbi',
    params:
        bed_prefix = 'results/eqtl/{datasource}/{group}/qqnorm.sorted',
    shell:
        '''
        bedtools sort -header -i {input} | bgzip -c > {output.bed}
        tabix -p bed {output.bed}

        # split file by chromosome
        for i in {{1..22}}; do
            tabix --print-header {output.bed} chr$i| bgzip -c > {params.bed_prefix}.chr$i.bed.gz
            tabix -p bed {params.bed_prefix}.chr$i.bed.gz
        done

        '''

rule GTExPhenoPC:
    input: 
        bed = 'results/eqtl/{datasource}/{group}/qqnorm.sorted.{chrom}.bed.gz',
    output: 'results/eqtl/{datasource}/{group}/qqnorm.sorted.{chrom}.pca'
    params:
        R_script = 'workflow/scripts/fromBen_PermuteAndPCA.R',
    shell:
        '''
        Rscript {params.R_script} {input} {output}
        '''

rule MakeCovMatrixEqtl:
    input:
        samples = 'results/pheno/noisy/{datasource}/{group}/separateNoise/leafcutter_names.txt',  # for changes
        GenoPCs = 'results/geno/{datasource}/{group}/{chrom}.pca',
        PhenoPCs = 'results/eqtl/{datasource}/{group}/qqnorm.sorted.{chrom}.pca',
    output: 'results/eqtl/{datasource}/{group}/qqnorm.sorted.{chrom}.cov'
    params:
        n_phenoPCs = 11, # index starts from header row
        n_genoPCs = 5, # index starts from header row
    run:
        fout = open(output[0], 'w')
        with open(input.samples) as f:
            samples = f.readlines()
            samples = [x.strip() for x in samples]
            fout.write('\t'.join(['id'] + samples) + '\n')

        with open(input.GenoPCs) as f: # append genotype PCs
            genoHeader = f.readline().strip().split()
            print("Getting genotype PCs...")
            if not all(a == b for a,b in zip(samples, genoHeader[1:])):
                print(f"Samples in genotype do not match!\nsamples:{samples}\nGenotype:{genoHeader[1:]}")
                exit("Samples genotype do not match! Exiting...")
            i = 1
            for ln in f.readlines():
                if i > params.n_genoPCs:
                    break
                ln = ln.strip().split()
                ln[0] = f'genoPC:{ln[0]}'
                fout.write('\t'.join(ln) + '\n')
                i += 1

        with open(input.PhenoPCs) as f: # append phenotype PCs
            print("Getting phenotype PCs...")
            phenoHeader = f.readline().strip().split()
            if not all(a == b for a,b in zip(samples, phenoHeader[1:])):
                exit("Samples in phenotype do not match! Exiting...")
            i = 1
            for ln in f.readlines():
                if i > params.n_phenoPCs:
                    break
                ln = ln.strip().split()
                ln[0] = f'phenoPC:{ln[0]}'
                fout.write('\t'.join(ln) + '\n')
                i += 1

        print(f"Done.\nWrote to covariance matrix file: {output[0]}")
        fout.close()




# ----------------------------------------------------------------------------------------
# run permutation pass eQTL for GTEx
# ----------------------------------------------------------------------------------------

rule GTEx_eQTL_Perm:
    message: 'Run permutation pass eQTL for GTEx'
    input: 
        vcf = 'results/geno/{datasource}/{group}/{chrom}.vcf.gz',
        cov = 'results/eqtl/{datasource}/{group}/qqnorm.sorted.{chrom}.cov',
    output: temp('results/eqtl/{datasource}/{group}/perm/{chrom}.txt')
    params:
        bed = 'results/eqtl/{datasource}/{group}/qqnorm.sorted.{chrom}.bed.gz',
        window = 1000000, # 1Mb
        chrom = '{chrom}'
    # log: 'logs/GTEx_eQTL_Perm_{datasource}_{group}_{chrom}.log'
    shell:
        '''
        module unload gsl && module load gsl/2.5

        QTLtools cis \
            --seed 123 \
            --permute 1000 \
            --vcf {input.vcf} --bed {params.bed} --cov {input.cov}  --out {output} \
            --window {params.window}
        '''

use rule AddQvalueToPermutationPass as AddQvalueToGTExEQtlPerm with:
    input:  'results/eqtl/{datasource}/{group}/perm/{chrom}.txt'
    output: 'results/eqtl/{datasource}/{group}/perm/{chrom}.txt.addQval.txt.gz'
    log: 'logs/AddQvalueToGTExEQtlPerm_{datasource}_{group}_{chrom}.log'


# ----------------------------------------------------------------------------------------
# run nominal pass eQTL for GTEx - simple (not suitable for coloc)
# simple version is meant for plotting sqtl-eqtl beta-beta correlation
# it does not include all snps, and didn't include genotype covariants
# ----------------------------------------------------------------------------------------

# Run nominal pass eQTL on top sQTL variants, used for plotting sQTL-eQTL beta-beta correlation
rule GTEx_eQTL_Nominal:
    message: 'Run nominal pass eQTL for GTEx, with only sQTLs top variants'
    input: 
        vcf = 'results/qtl/noisy/{datasource}/{group}/separateNoise/cis_100000/nom/{chrom}.TopVariants.vcf.gz', # top varaints from sQTL
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


# Run nominal pass eQTL for GTEx with all genomewide SNPs (cis), not just top sQTL variants
# Use it for colocalization analysis
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




