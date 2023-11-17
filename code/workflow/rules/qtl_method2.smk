# NOTE qtltools require gsl, make sure gsl is loaded before running snakemake

# -------------------- All QTL related rules --------------------



def getExtractNoisyPhenoInput(wildcards):
    if wildcards.datasource == 'GTEx':
        counts = os.path.join('results/pheno/noisy', wildcards.datasource, 
                wildcards.group, 'leafcutter_perind.counts.noise.gz')
        vcf = config['annotation']['GTEx']['VCF_v7_HG38']
    elif wildcards.datasource == 'Geuvadis':
        counts = os.path.join('results/pheno/noisy', wildcards.datasource, 
                wildcards.group, 'leafcutter_perind.counts.noise.gz')
        vcf = '/project2/yangili1/zpmu/1kg_b38/CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz'
        # only need the vcf header to get all available sample names, thus choose a smaller vcf
    return {'counts': counts, 'vcf': vcf}

rule ExtractNoisyCounts:
    input: unpack(getExtractNoisyPhenoInput)
    output:
        counts = 'results/pheno/noisy/{datasource}/{group}/leafcutter_perind.counts.noiseonly.gz',
        indivs = 'results/pheno/noisy/{datasource}/{group}/individuals.txt'
    run:
        import subprocess
        vcf_header = f"bcftools view -h {input.vcf} | awk '$1 ~ /#CHROM/' "
        vcf_header = subprocess.run(vcf_header, shell=True, capture_output=True).stdout.decode('utf-8').strip().split()
        geno_inds = vcf_header[9:] # available genotype individual IDs

        with gzip.open(input.counts) as f:
            fout = gzip.open(output.counts, 'wt')
            i = 0
            for ln in f:
                lnsplit = ln.decode().strip().split()
                if i == 0: # header
                    keep1 = [lnsplit.index(l) for l in lnsplit if l in ['chrom'] + geno_inds ] # cols to keep
                    keep2 = [lnsplit.index(l) for l in lnsplit if l in geno_inds ] # cols to keep, no 'chrom'
                    with open(output.indivs, 'w') as sf: # write individual ids
                        indvs = [lnsplit[ix] for ix in keep2]
                        print(f'Removing these samples not in genotype vcf:\n{[x for x in lnsplit[1:] if x not in indvs]} ...')
                        sf.writelines([x + '\n' for x in indvs])
                    buf = ' '.join([lnsplit[ix] for ix in keep1]) + '\n' # header line
                    fout.write(buf)
                elif '*' in lnsplit[0]:
                    buf = ' '.join([lnsplit[ix] for ix in keep1]) + '\n'
                    fout.write(buf)
                i += 1


rule PrepPhenoTable:
    input: rules.ExtractNoisyCounts.output.counts
    output: 'results/pheno/noisy/{datasource}/{group}/leafcutter_perind.counts.noiseonly.gz.PCs'
    params: 
        py_script = 'workflow/submodules/leafcutter/scripts/prepare_phenotype_table.py',
        nPCs = 5
    conda: 'leafcutter'
    shell: 
        '''

        python {params.py_script} -p {params.nPCs} {input}
        ls {output}
        '''



rule MakePhenoBed:
    message:'### Make QTLtools required phenotype bed file for {wildcards.datasource}:{wildcards.group}:{wildcards.chrom}'
    input:
        PCs = rules.PrepPhenoTable.output
    output:
        bed = 'results/pheno/noisy/{datasource}/{group}/method2/{chrom}.bed.gz'
    params:
        out_dir = 'results/pheno/noisy/{datasource}/{group}/method2',
        pheno = 'results/pheno/noisy/{datasource}/{group}/leafcutter_perind.counts.noiseonly.gz.phen_{chrom}'
    run:
        outname = f'{output.bed}'.replace('.bed.gz', '')
        with open(outname + '.bed', 'w') as fout:
            with open(params.pheno) as f:
                i = 0
                for ln in f:
                    lnsplit = ln.strip().split()
                    datacols = lnsplit[4:]
                
                    if i == 0:
                        idcols = ['#Chr', 'start', 'end', 'pid', 'gid', 'strand']
                    else:
                        chrom, start, end, pid = lnsplit[0:4]
                        gid = pid
                        strand = pid.split("_")[2]
                        idcols = [chrom, start, end, pid, gid, strand]
                    
                    buf = '\t'.join(idcols + datacols) + '\n'
                    fout.write(buf)
                    
                    i += 1
        print(f'Wrote {i} lines ...')    
        bgzip = f'cat <(head -1 {outname}.bed) <(sortBed -i {outname}.bed) | bgzip -c > {output.bed}'
        print(f'run {bgzip} ...')
        shell(bgzip)
        tabix = f'tabix -p bed {output.bed}'
        print(f'run {tabix} ...')
        shell(tabix)
        print('Done.')



def getExtractGenotypeInput(wildcards):
    if wildcards.datasource == 'GTEx':
        vcf = config['annotation']['GTEx']['VCF_v7_HG38']
    elif wildcards.datasource == 'Geuvadis':
        vcf =  f'/project2/yangili1/zpmu/1kg_b38/CCDG_14151_B01_GRM_WGS_2020-08-05_{wildcards.chrom}.filtered.shapeit2-duohmm-phased.vcf.gz',
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
        PhenoPCs = 'results/pheno/noisy/{datasource}/{group}/leafcutter_perind.counts.noiseonly.gz.PCs',
        GenoPCs = rules.GenotypePCA.output
    output: 'results/pheno/noisy/{datasource}/{group}/method2/{chrom}_CovMatrix.txt'
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
        vcf = 'results/geno/{datasource}/{group}/{chrom}.vcf.gz',
        bed = 'results/pheno/noisy/{datasource}/{group}/method2/{chrom}.bed.gz',
        cov = 'results/pheno/noisy/{datasource}/{group}/method2/{chrom}_CovMatrix.txt'
    output: temp('results/qtl/noisy/{datasource}/{group}/cis_{window}/perm/method2/{chrom}.txt')
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


rule AddQvalueToPermutationPass:
    message: '### Add qvalue to the output of QTLtools permutation pass'
    input:  'results/qtl/noisy/{datasource}/{group}/cis_{window}/perm/method2/{chrom}.txt'
    output: 'results/qtl/noisy/{datasource}/{group}/cis_{window}/perm/method2/{chrom}.addQval.txt.gz'
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






# this rule is obsolete
# rule ExtractPhenotypeBySamples:
#     '''
#     Procedures:
#         1.  convert ERR ID in column into Sample ID
#         2.  select only columns linked with genotype, and having 1 ERR per sample
#         3.  select only rows representing noisy splicing (with *)
#         4.  output noisy count table with selected rows and columns
#         5.  output selected sample ID in a file for genotype prep procedure
#     '''
#     message: 'Prepare phenotype count table'
#     input: 
#         counts = 'results/pheno/Geuvadis/noise_counts/Geuvadis_perind.counts.noise.gz',
#         metadata = config['Dataset']['Geuvadis']['Metadata'],
#         Linked_1to1_Samples = config['Dataset']['Geuvadis']['Linked_1to1_SampleIDs']
#     output: 
#         counts     = 'results/pheno/{datasource}/{population}/noise_pheno/pheno_extract.count.gz',
#         samplelist = 'results/pheno/{datasource}/{population}/noise_pheno/samplelist_pheno.txt'
#     params: 
#         py_script = 'workflow/scripts/extract_pheno.py',
#         population = "{population}"
#     threads: 1
#     resources: cpu = 1, mem_mb = 15000, time = 2100
#     shell: 
#         '''
#             python {params.py_script} -I {input.counts} \
#                 --outcount {output.counts} \
#                 --outsample {output.samplelist} \
#                 --lookuptable {input.metadata} \
#                 --subset {input.Linked_1to1_Samples} \
#                 --pop {params.population}
#         '''
    # in future analysis, use geuvadis-1kgp-common-sample-id.txt to subset samples
    # such that all Geuvadis samples with a 1KGP genotype match is selected.
    # Currently, only samples with only 1 ERR id and with a match with 1kgp genotype is selected.
    # NOTE by this subselect, YRI only has 21 samples with genotype and having 1 ERR_id only.
    # thus for this practice run, only run EUR sanokes,




