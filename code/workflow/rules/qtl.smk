'''
This script is intended to work with both GTEx datasets and Geuvadis datasets.

datasource: GTEx or Geuvadis
group: tissue for GTEx, e.g. Bladder, or population for Geuvadis, e.g. EUR

'''

def getPreparePhenoBedInput(wildcards):
    datasource, group, phenType = wildcards.datasource, wildcards.group, wildcards.phenType
    if phenType == 'separateNoise':
        return 'results/pheno/noisy/{datasource}/{group}/leafcutter_perind.counts.noise_by_intron.gz'
    elif phenType == 'combineNoise':
        return 'results/pheno/noisy/{datasource}/{group}/leafcutter_perind.counts.noise.gz'
    else:
        print('Error. Invalid phenotype type. Exiting...')
        exit(1)

def getVcfIndiv(wildcards):
    '''get individual ids of the vcf
    '''
    datasource, group = wildcards.datasource, wildcards.group
    if datasource == 'GTEx':
        return config['VCF']['GTEx']['HG38_v7_indivs']
    elif datasource == 'Geuvadis':
        return config['VCF']['Geuvadis']['HG38_1kg_b38_indivs']
    else:
        print('Error. Invalid datasource. Exiting...')
        exit(1)

rule PreparePhenoBed:
    message: '### Prepare phenotype bed file for qtltools, using *.counts.noise.gz'
    input: getPreparePhenoBedInput
    output: 
        flag = touch('results/pheno/noisy/{datasource}/{group}/{phenType}/done'),
        samples = 'results/pheno/noisy/{datasource}/{group}/{phenType}/leafcutter_names.txt'
    params:
        pyscript = 'workflow/scripts/preparePheno.py',
        outPrefix = 'results/pheno/noisy/{datasource}/{group}/{phenType}/leafcutter',
        vcfSamples = getVcfIndiv ,# individual ids in vcf file
    log: 'logs/PreparePhenoBed_{datasource}_{group}_{phenType}.log'
    shell:
        '''
        python {params.pyscript} {input} --sampleFile {params.vcfSamples} --outPrefix {params.outPrefix} &> {log}

        for i in {{1..22}}; do
            bgzip -f {params.outPrefix}.phen_chr${{i}} 2>> {log}
            bgzip -f {params.outPrefix}.qqnorm_chr${{i}} 2>> {log}
            tabix -f -p bed {params.outPrefix}.qqnorm_chr${{i}}.gz 2>> {log}
        done

        ls -l {output.samples} &>> {log}


        '''


def getExtractGenotypeInput(wildcards):
    if wildcards.datasource == 'GTEx':
        vcf = config['VCF']['GTEx']['HG38_v7']
    elif wildcards.datasource == 'Geuvadis':
        vcf_dir = config['VCF']['Geuvadis']['HG38_1kg_b38']
        chrom = wildcards.chrom
        vcf =  f'{vcf_dir}/CCDG_14151_B01_GRM_WGS_2020-08-05_{chrom}.filtered.shapeit2-duohmm-phased.vcf.gz'
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
    message: '### Extract genotype vcf for phenotype samples'
    input: 
        vcf = getExtractGenotypeInput,
        sample_file = 'results/pheno/noisy/{datasource}/{group}/separateNoise/leafcutter_names.txt'
    output: 
        vcf = 'results/geno/{datasource}/{group}/{chrom}.vcf.gz',
        tbi = 'results/geno/{datasource}/{group}/{chrom}.vcf.gz.tbi'
    params:
        expr = getExtractGenotypeParams
    log: 'logs/ExtractGenotypeVCF_{datasource}_{group}_{chrom}.log'
    threads: 4
    resources: time=2000, mem_mb=15000, cpu=4
    group: 'geno'
    shell:
        '''
        bcftools view \
            --threads {threads} \
            --samples-file {input.sample_file} \
            -i "{params.expr}" \
            -Oz -o {output.vcf} \
            {input.vcf} \
            {wildcards.chrom} &> {log}
        
        bcftools index --threads {threads} --tbi {output.vcf} &>> {log}
        '''

# sometimes qtltools fail, module unload gsl gcc and reload gcc/10.2.0 and gsl/2.2.1 works
rule GenotypePCA:
    input:  'results/geno/{datasource}/{group}/{chrom}.vcf.gz'
    output: 'results/geno/{datasource}/{group}/{chrom}.pca'
    params:
        out_prefix = 'results/geno/{datasource}/{group}/{chrom}'
    log: 'logs/GenotypePCA_{datasource}_{group}_{chrom}.log'
    group: 'geno'
    shell:
        '''
            module unload gsl && module load gsl/2.5
            QTLtools pca \
            --seed 123 \
            --maf 0.05 \
            --vcf {input}  \
            --out {params.out_prefix} \
            --center \
            --scale &> {log}
        '''


rule MakeCovarianceMatrix:
    message: '### Make covariance matrix with fixed PCs'
    input:
        samples = 'results/pheno/noisy/{datasource}/{group}/{phenType}/leafcutter_names.txt',  # for changes
        GenoPCs = 'results/geno/{datasource}/{group}/{chrom}.pca',
    output: 'results/pheno/noisy/{datasource}/{group}/{phenType}/{chrom}_CovMatrix.txt'
    params:
        n_phenoPCs = 11, # index starts from header row
        n_genoPCs = 5, # index starts from header row
        PhenoPCs = 'results/pheno/noisy/{datasource}/{group}/{phenType}/leafcutter.PCs',
        fake = "fake",
    log: 'logs/MakeCovarianceMatrix_{datasource}_{group}_{phenType}_{chrom}.log'
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

        with open(params.PhenoPCs) as f: # append phenotype PCs
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



rule MapQTL_Perm:
    '''
        already ranknorm normalized, do not need the normal option 
    '''
    message: 'Map QTL using permutation pass'
    input:
        phenoPrep = 'results/pheno/noisy/{datasource}/{group}/{phenType}/done',
        vcf = 'results/geno/{datasource}/{group}/{chrom}.vcf.gz',
        cov = 'results/pheno/noisy/{datasource}/{group}/{phenType}/{chrom}_CovMatrix.txt',
    output: temp('results/qtl/noisy/{datasource}/{group}/{phenType}/cis_{window}/perm/{chrom}.txt')
    log: 'logs/MapQTL_Perm_{datasource}_{group}_{phenType}_{window}_{chrom}.log'
    params:
        cis_window = '{window}',
        pheno = 'results/pheno/noisy/{datasource}/{group}/{phenType}/leafcutter.qqnorm_{chrom}.gz',
    resources: cpu = 1, mem = 12000, time = 1000
    shell:
        '''
        module unload gsl && module load gsl/2.5
        QTLtools cis \
            --seed 123 \
            --vcf {input.vcf} --bed {params.pheno} --cov {input.cov}  --out {output} \
            --window {params.cis_window} \
            --permute 1000 &> {log}
        '''


rule AddQvalueToPermutationPass:
    message: '### Add qvalue to the output of QTLtools permutation pass'
    input:  'results/qtl/noisy/{datasource}/{group}/{phenType}/cis_{window}/perm/{chrom}.txt'
    output: 'results/qtl/noisy/{datasource}/{group}/{phenType}/cis_{window}/perm/{chrom}.addQval.txt.gz'
    params:
        rscript = 'workflow/scripts/AddQvalueToQTLtoolsOutput.R'
    shell:
        '''
            txt=$(echo {output} | sed -E 's/.gz//')
            Rscript {params.rscript} {input} $txt
            bgzip -f $txt
        '''
#     # QTLtools cis permutation pass output fields:
#     # - 1. phenotype_id
#     # - 2. phenotype_chr
#     # - 3. phenotype_start
#     # - 4. phenotype_end
#     # - 5. phenotype_strand
#     # - 6. number_variants_tested
#     # - 7. best_nominal_distance
#     # - 8. best_genotype_id
#     # - 9. best_genotype_chr
#     # - 10. best_genotype_start
#     # - 11. best_genotype_end
#     # - 12. true_degrees_of_freedom
#     # - 13. estimated_dof
#     # - 14. beta_ml1
#     # - 15. beta_ml2
#     # - 16. pval_nominal
#     # - 17. pval_r2
#     # - 18. regression_slop
#     # - 19. pval_empirical (by direct permutation method)
#     # - 20. pval_adjusted (adjusted p-value, note, this has not been genome-wide multiple-tested on phenotypes)


rule ExtractSNPsForNominal:
    message: '### Extract top SNPs of sQTL for nominal pass'
    input: 
        perm = 'results/qtl/noisy/{datasource}/{group}/{phenType}/cis_{window}/perm/{chrom}.addQval.txt.gz',
        vcf = 'results/geno/{datasource}/{group}/{chrom}.vcf.gz'
    output: 
        bed = 'results/qtl/noisy/{datasource}/{group}/{phenType}/cis_{window}/nom/{chrom}.TopVariants.bed',
        vcf = 'results/qtl/noisy/{datasource}/{group}/{phenType}/cis_{window}/nom/{chrom}.TopVariants.vcf.gz'
    log: 'logs/ExtractSNPsForNominal_{datasource}_{group}_{phenType}_{window}_{chrom}.log'
    shell:
        '''
        (zcat {input.perm} | \
            awk 'BEGIN {{OFS="\t"}};
                 NR>1 {{print $9, $10 - 1, $11}}
                ' | \
            sort -k2 -k3 | uniq) 1> {output.bed} 2> {log}

        bcftools view \
            -R {output.bed} \
            -Oz -o {output.vcf} \
            {input.vcf} &>> {log}
        
        bcftools index --tbi {output.vcf} &>> {log}

        '''

rule NominalQTL:
    message: '### Map QTL using nominal pass'
    input: 
        phenoPrep = 'results/pheno/noisy/{datasource}/{group}/{phenType}/done',
        vcf = 'results/qtl/noisy/{datasource}/{group}/{phenType}/cis_{window}/nom/{chrom}.TopVariants.vcf.gz',
        cov = 'results/pheno/noisy/{datasource}/{group}/{phenType}/{chrom}_CovMatrix.txt',
    output: temp('results/qtl/noisy/{datasource}/{group}/{phenType}/cis_{window}/nom/{chrom}.txt')
    log: 'logs/NominalQTL_{datasource}_{group}_{phenType}_{window}_{chrom}.log'
    params:
        cis_window = lambda w: int(w.window) + 10000, # add 10kb to the window
        pheno = 'results/pheno/noisy/{datasource}/{group}/{phenType}/leafcutter.qqnorm_{chrom}.gz',
    resources: cpu = 1, mem = 12000, time = 1000
    shell:
        '''
        module unload gsl && module load gsl/2.5
        QTLtools cis \
            --seed 123 \
            --nominal 1 \
            --vcf {input.vcf} --bed {params.pheno} --cov {input.cov}  --out {output} \
            --window {params.cis_window} &> {log}
        '''

rule TabixNominal:
    message: 'Tabix the nominal pass output'
    input: 'results/qtl/noisy/{datasource}/{group}/{phenType}/cis_{window}/nom/{chrom}.txt'
    output: 'results/qtl/noisy/{datasource}/{group}/{phenType}/cis_{window}/nom/{chrom}.txt.gz'
    log: 'logs/TabixNominal_{datasource}_{group}_{phenType}_{window}_{chrom}.log'
    shell:
        '''
        (awk 'BEGIN {{OFS="\t"}}; 
                    {{for (i=1; i<=NF; i++) printf "%s%s", $i, (i<NF ? OFS : ORS)}}
             ' {input} | \
            sort -k9 -k10 -V | \
            bgzip -c > {output}) 2> {log}
        (tabix -s 9 -b 10 -e 11 {output}) &> {log}
        '''

rule sortGTExPhenotype_rawPSI:
    input: 'results/pheno/noisy/{datasource}/{group}/{phenType}/leafcutter.phen_{chrom}.gz',
    output: 
        bed = 'results/pheno/noisy/{datasource}/{group}/{phenType}/leafcutter.phen_{chrom}_sorted.gz',
        tbi = 'results/pheno/noisy/{datasource}/{group}/{phenType}/leafcutter.phen_{chrom}_sorted.gz.tbi',
    shell:
        '''
        bedtools sort -header -i {input} | bgzip -c > {output.bed}
        tabix -p bed {output.bed}

        '''


rule GTExRawPsiPca:
    message: """run PCA on raw PSI"""
    input: 
        bed = 'results/pheno/noisy/{datasource}/{group}/{phenType}/leafcutter.phen_{chrom}_sorted.gz',
    output:
        pca = 'results/pheno/noisy/{datasource}/{group}/{phenType}/leafcutter.phen_{chrom}.pca',
    params:
        R_script = 'workflow/scripts/fromBen_PermuteAndPCA.R'
    shell:
        '''
        Rscript {params.R_script} {input} {output}
        '''


rule GTEx_sQTL_Nominal_rawPSI:
    message: """Run nominal pass sQTL for GTEx, with only sQTLs top variants and raw psi values"""
    input: 
        vcf = 'results/qtl/noisy/{datasource}/{group}/separateNoise/cis_100000/nom/{chrom}.TopVariants.vcf.gz',
        bed = 'results/pheno/noisy/{datasource}/{group}/separateNoise/leafcutter.phen_{chrom}_sorted.gz',
        cov = 'results/pheno/noisy/{datasource}/{group}/separateNoise/leafcutter.phen_{chrom}.pca',
    output: 'results/qtl/noisy/{datasource}/{group}/separateNoise/cis_100000/nom-raw-psi/{chrom}.txt'
    params:
        window = 1000000, # 1Mb, doesn't matter, since vcf has only top variants
        chrom = '{chrom}'
    log: 'logs/GTEx_sQTL_Nominal_rawPSI_{datasource}_{group}_{chrom}.log'
    shell:
        '''
        module unload gsl && module load gsl/2.5
        QTLtools cis \
            --seed 123 \
            --nominal 1 \
            --vcf {input.vcf} --bed {input.bed} --cov {input.cov}  --out {output} \
            --window {params.window} &> {log}
        '''



#-----------------------------------------------------------------
#            Test pseudo count method
#-----------------------------------------------------------------


# test qtl tools with different pseudocount method
use rule PreparePhenoBed as PreparePhenoBed_test with:
    message: '### Prepare phenotype bed file for qtltools, different pseudocount method'
    input: 'results/pheno/noisy/{datasource}/{group}/leafcutter_perind.counts.noise_by_intron.gz'
    output:
        flag = touch('results/pheno/noisy/{datasource}/{group}/test_pseudocount/done'),
        samples = 'results/pheno/noisy/{datasource}/{group}/test_pseudocount/leafcutter_names.txt'
    params:
        pyscript = 'workflow/scripts/preparePheno-test-pseudocount.py',
        outPrefix = 'results/pheno/noisy/{datasource}/{group}/test_pseudocount/leafcutter',
        vcfSamples = getVcfIndiv ,# individual ids in vcf file
    log: 'logs/PreparePhenoBed_{datasource}_{group}_test_pseudocount.log'


use rule MakeCovarianceMatrix as MakeCovarianceMatrix_test with:
    message: '### Make covariance matrix with fixed PCs, different pseudocount method'
    input:
        samples = 'results/pheno/noisy/{datasource}/{group}/test_pseudocount/leafcutter_names.txt', 
        GenoPCs = 'results/geno/{datasource}/{group}/{chrom}.pca',
        flag = 'results/pheno/noisy/{datasource}/{group}/test_pseudocount/done',
    output: 'results/pheno/noisy/{datasource}/{group}/test_pseudocount/{chrom}_CovMatrix.txt'
    params:
        n_phenoPCs = 11, # index starts from header row
        n_genoPCs = 5, # index starts from header row
        PhenoPCs = 'results/pheno/noisy/{datasource}/{group}/test_pseudocount/leafcutter.PCs',
        fake = "fake",
    log: 'logs/MakeCovarianceMatrix_{datasource}_{group}_test_pseudocount_{chrom}.log'


use rule MapQTL_Perm as MapQTL_Perm_test with:
    message: 'Map QTL using permutation pass, different pseudocount method'
    input:
        phenoPrep = 'results/pheno/noisy/{datasource}/{group}/test_pseudocount/done',
        vcf = 'results/geno/{datasource}/{group}/{chrom}.vcf.gz',
        cov = 'results/pheno/noisy/{datasource}/{group}/test_pseudocount/{chrom}_CovMatrix.txt',
    output: temp('results/qtl/noisy/{datasource}/{group}/test_pseudocount/cis_{window}/perm/{chrom}.txt')
    log: 'logs/MapQTL_Perm_{datasource}_{group}_test_pseudocount_{window}_{chrom}.log'
    params:
        cis_window = '{window}',
        pheno = 'results/pheno/noisy/{datasource}/{group}/test_pseudocount/leafcutter.qqnorm_{chrom}.gz',
    resources: cpu = 1, mem = 12000, time = 1000


use rule AddQvalueToPermutationPass as AddQvalueToPermutationPass_test with:
    message: '### Add qvalue to the output of QTLtools permutation pass, different pseudocount method'
    input:  'results/qtl/noisy/{datasource}/{group}/test_pseudocount/cis_{window}/perm/{chrom}.txt'
    output: 'results/qtl/noisy/{datasource}/{group}/test_pseudocount/cis_{window}/perm/{chrom}.addQval.txt.gz'
    log: 'logs/AddQvalueToPermutationPass_{datasource}_{group}_test_pseudocount_{window}_{chrom}.log'


use rule ExtractSNPsForNominal as ExtractSNPsForNominal_test with:
    message: '### Extract top SNPs of sQTL for nominal pass, different pseudocount method'
    input: 
        perm = 'results/qtl/noisy/{datasource}/{group}/test_pseudocount/cis_{window}/perm/{chrom}.addQval.txt.gz',
        vcf = 'results/geno/{datasource}/{group}/{chrom}.vcf.gz'
    output:
        bed = 'results/qtl/noisy/{datasource}/{group}/test_pseudocount/cis_{window}/nom/{chrom}.TopVariants.bed',
        vcf = 'results/qtl/noisy/{datasource}/{group}/test_pseudocount/cis_{window}/nom/{chrom}.TopVariants.vcf.gz'
    log: 'logs/ExtractSNPsForNominal_{datasource}_{group}_test_pseudocount_{window}_{chrom}.log'


use rule NominalQTL as NominalQTL_test with:
    message: '### Map QTL using nominal pass, different pseudocount method'
    input: 
        phenoPrep = 'results/pheno/noisy/{datasource}/{group}/test_pseudocount/done',
        vcf = 'results/qtl/noisy/{datasource}/{group}/test_pseudocount/cis_{window}/nom/{chrom}.TopVariants.vcf.gz',
        cov = 'results/pheno/noisy/{datasource}/{group}/separateNoise/{chrom}_CovMatrix.txt',
    output: temp('results/qtl/noisy/{datasource}/{group}/test_pseudocount/cis_{window}/nom/{chrom}.txt')
    log: 'logs/NominalQTL_{datasource}_{group}_test_pseudocount_{window}_{chrom}.log'
    params:
        cis_window = lambda w: int(w.window) + 10000, # add 10kb to the window
        pheno = 'results/pheno/noisy/{datasource}/{group}/test_pseudocount/leafcutter.qqnorm_{chrom}.gz',
    resources: cpu = 1, mem = 12000, time = 1000

use rule TabixNominal as TabixNominal_test with:
    message: 'Tabix the nominal pass output, different pseudocount method'
    input: 'results/qtl/noisy/{datasource}/{group}/test_pseudocount/cis_{window}/nom/{chrom}.txt'
    output: 'results/qtl/noisy/{datasource}/{group}/test_pseudocount/cis_{window}/nom/{chrom}.txt.gz'
    log: 'logs/TabixNominal_{datasource}_{group}_test_pseudocount_{window}_{chrom}.log'




