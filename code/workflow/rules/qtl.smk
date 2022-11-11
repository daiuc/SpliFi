# -------------------- All QTL related rules --------------------


rule ExtractPhenotypeBySamples:
    '''
    Procedures:
        1.  convert ERR ID in column into Sample ID
        2.  select only columns linked with genotype, and having 1 ERR per sample
        3.  select only rows representing noisy splicing (with *)
        4.  output noisy count table with selected rows and columns
        5.  output selected sample ID in a file for genotype prep procedure
    '''
    message: 'Prepare phenotype count table'
    input: 
        counts = 'results/pheno/Geuvadis/noise_counts/Geuvadis_perind.counts.noise.gz',
        metadata = config['Dataset']['Geuvadis']['Metadata'],
        Linked_1to1_Samples = config['Dataset']['Geuvadis']['Linked_1to1_SampleIDs']
    output: 
        counts     = 'results/pheno/{datasource}/{population}/noise_pheno/pheno_extract.count.gz',
        samplelist = 'results/pheno/{datasource}/{population}/noise_pheno/samplelist_pheno.txt'
    params: 
        py_script = 'workflow/scripts/extract_pheno.py',
        population = "{population}"
    threads: 1
    resources: cpu = 1, mem_mb = 15000, time = 2100
    shell: 
        '''
            python {params.py_script} -I {input.counts} \
                --outcount {output.counts} \
                --outsample {output.samplelist} \
                --lookuptable {input.metadata} \
                --subset {input.Linked_1to1_Samples} \
                --pop {params.population}
        '''
    # in future analysis, use geuvadis-1kgp-common-sample-id.txt to subset samples
    # such that all Geuvadis samples with a 1KGP genotype match is selected.
    # Currently, only samples with only 1 ERR id and with a match with 1kgp genotype is selected.
    # NOTE by this subselect, YRI only has 21 samples with genotype and having 1 ERR_id only.
    # thus for this practice run, only run EUR sanokes,


rule ExtractGenotypeVCF:
    '''
        Procedures:
            1.  extract vcf per chromosome for samples in phenotype
            2.  tabix vcf
    '''
    message: '### Extract genotype vcf for phenotype samples.'
    input: 
        vcf =  '/project2/yangili1/zpmu/1kg_b38/CCDG_14151_B01_GRM_WGS_2020-08-05_{chrom}.filtered.shapeit2-duohmm-phased.vcf.gz',
        sample_file = rules.ExtractPhenotypeBySamples.output.samplelist
    output: 'results/geno/{datasource}/{population}/{chrom}.vcf.gz'
    params:
        min_MAF = 0.05,
        max_HWE = 1e-3
    threads: 4
    resources: time=2000, mem_mb=15000, cpu=4
    wildcard_constraints:
        chrom = "chr[0-9]{1,2}"
    shell:
        '''
        bcftools view \
            --threads {threads} \
            --samples-file {input.sample_file} \
            -i "AF >= {params.min_MAF} && HWE < {params.max_HWE}" \
            -Oz -o {output} \
            {input.vcf}
        
        bcftools index --threads {threads} --tbi {output}
        '''

rule PrepPhenoTable:
	input: rules.ExtractPhenotypeBySamples.output.counts
	output: touch('results/pheno/{datasource}/{population}/noise_pheno/PrepPhenoTable.done')
	params:
		py_script = 'workflow/submodules/leafcutter/scripts/prepare_phenotype_table.py'
	conda: 'leafcutter'
	shell:
		'''
		python {params.py_script} {input}
		'''

rule MakePhenoBed:
    message:'### Make phenotype bed format required by QTLtools'
    input:
        flag  = 'results/pheno/{datasource}/{population}/noise_pheno/PrepPhenoTable.done',
        pheno = 'results/pheno/{datasource}/{population}/noise_pheno/pheno_extract.count.gz.phen_{chrom}'
    output:
        pheno = 'results/pheno/{datasource}/{population}/noise_pheno/{chrom}.bed.gz'
    params:
        rscript = 'workflow/scripts/Make_pheno_bed.R',
        out_dir = 'results/pheno/{datasource}/{population}/noise_pheno/'
    threads: 1
    shell:
        '''
        temp_out={params.out_dir}$(basename -s .gz {output.pheno})
        Rscript {params.rscript} {input.pheno} $temp_out
        bgzip $temp_out
        tabix -p bed {output.pheno}
        '''


rule PhenotypePCA:
    message: '### Run PCA on phenotype with permutation'
    input: rules.MakePhenoBed.output.pheno
    output: 'results/pheno/{datasource}/{population}/noise_pheno/{chrom}.pca'
    params:
        rscript = 'workflow/scripts/PermuteAndPCA.R'
    threads: 1
    shell:
        '''
        Rscript {params.rscript} {input} {output} 
        '''

rule GenotypePCA:
    message: '### Run pca on genotype'
    input:  'results/geno/{datasource}/{population}/{chrom}.vcf.gz'
    output: 'results/geno/{datasource}/{population}/{chrom}.pca'
    params:
        out_prefix = 'results/geno/{datasource}/{population}/{chrom}'
    threads: 1
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
    output: 'results/pheno/{datasource}/{population}/noise_pheno/{chrom}_CovMatrix.txt'
    params:
        Geno_PCs = 4,
        rscript = 'workflow/scripts/Make_CovMatrix.R'
    shell:
        '''
        Rscript {params.rscript} {input.PhenoPCs} {input.GenoPCs} {params.Geno_PCs} {output}
        '''


rule MapQTL_Perm:
    message: 'Map QTL using permutation pass'
    input:
        vcf = 'results/geno/{datasource}/{population}/{chrom}.vcf.gz',
        bed = 'results/pheno/{datasource}/{population}/noise_pheno/{chrom}.bed.gz',
        cov = 'results/pheno/{datasource}/{population}/noise_pheno/{chrom}_CovMatrix.txt'
    output: 'results/qtl/noise/{datasource}/{population}/cis_{window}/perm/{chrom}.txt'
    params:
        cis_window = '{window}'
    resources: cpu = 1, mem = 12000, time = 1000
    shell:
        '''
        # already quantile normalized, do not need the normal option 
        QTLtools cis \
            --seed 123 --silent \
            --vcf {input.vcf} --bed {input.bed} --cov {input.cov}  --out {output} \
            --window {params.cis_window} \
            --permute 1000 \
        '''


rule AddQvalueToPermutationPass:
    message: '### Add qvalue to the output of QTLtools permutation pass'
    input:  'results/qtl/noise/{datasource}/{population}/cis_{window}/perm/{chrom}.txt'
    output: 'results/qtl/noise/{datasource}/{population}/cis_{window}/perm/{chrom}.addQval.txt.gz'
    params:
        rscript = 'workflow/scripts/AddQvalueToQTLtoolsOutput.R'
    shell:
        '''
            Rscript {params.rscript} {input} {output}
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



