

rule RandomizePhenoBed:
    message: 'qtl QC - randomize phenotype bed file'
    input: 'results/pheno/noisy/{datasource}/{group}/pheno.{chrom}.bed.gz'
    output: 'results/qtlQC/noisy/{datasource}/{group}/pheno.{chrom}.bed.gz'
    threads: 1
    resources: cpu = 1, mem_mb = 15000, time = 2100
    run:
        import pandas as pd
        import random

        df = pd.read_csv(input[0], sep='\t')
        idcols = df.columns[0:6] # id columns
        dacols = df.columns[6:] # data columns
        
        # randomized data column rows
        rowids = random.sample(range(df.shape[0]), df.shape[0])

        # save randomized data in a new df
        df2 = pd.concat([df[idcols].reset_index(drop=True),
                df[dacols].iloc[rowids].reset_index(drop=True)], axis=1)
        
        # write new dataframe
        bedfile = output[0].replace('.gz', '')
        df2.to_csv(bedfile, sep="\t", header=True, index=False)

        # gzip
        bgzip = f'bgzip {bedfile}'
        shell(bgzip)

        # tabix
        tabix = f'tabix -p bed {output[0]}'
        shell(tabix)


rule RandomPhenotypePCA:
    message: '### Run PCA on randomized phenotype bed'
    input: rules.RandomizePhenoBed.output
    output: 'results/qtlQC/noisy/{datasource}/{group}/{chrom}.pca'
    params:
        rscript = 'workflow/scripts/PermuteAndPCA.R',
    shell:
        '''
        Rscript {params.rscript} {input} {output} 
        '''


rule RandomCovarianceMatrix:
    message: '### Make covariance matrix for randomized phenotypes'
    input:
        PhenoPCs = rules.RandomPhenotypePCA.output,
        GenoPCs = rules.GenotypePCA.output
    output: 'results/qtlQC/noisy/{datasource}/{group}/{chrom}_CovMatrix.txt'
    params:
        Geno_PCs = 4,
        rscript = 'workflow/scripts/Make_CovMatrix.R'
    shell:
        '''
        Rscript {params.rscript} {input.PhenoPCs} {input.GenoPCs} {params.Geno_PCs} {output}
        '''


rule RandomMapQTL_Perm:
    message: 'Map QTL using permutation pass'
    input:
        vcf = 'results/geno/{datasource}/{group}/{chrom}.vcf.gz',
        bed = 'results/qtlQC/noisy/{datasource}/{group}/pheno.{chrom}.bed.gz',
        cov = 'results/qtlQC/noisy/{datasource}/{group}/{chrom}_CovMatrix.txt'
    output: temp('results/qtlQC/perm/{datasource}/{group}/cis_{window}/{chrom}.txt')
    params:
        cis_window = '{window}'
    resources: cpu = 1, mem = 12000, time = 1000
    shell:
        '''
        # already ranknorm normalized, do not need the normal option 
        QTLtools cis \
            --seed 123 --silent \
            --vcf {input.vcf} --bed {input.bed} --cov {input.cov}  --out {output} \
            --window {params.cis_window} \
            --permute 1000 --region {wildcards.chrom}
        '''


rule RandomAddQvalueToPermutationPass:
    message: '### Add qvalue to the output of QTLtools permutation pass'
    input:  rules.RandomMapQTL_Perm.output
    output: 'results/qtlQC/perm/{datasource}/{group}/cis_{window}/{chrom}.addQval.txt.gz'
    params:
        rscript = 'workflow/scripts/AddQvalueToQTLtoolsOutput.R'
    shell:
        '''
            txt=$(echo {output} | sed -E 's/.gz//')
            Rscript {params.rscript} {input} $txt
            bgzip $txt
        '''





