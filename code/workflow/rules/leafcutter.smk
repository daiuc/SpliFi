'''
Use leafcutter to cluster junctions into splicing events
Then use qtltools to perform QTL mapping

'''

configfile: "config/config.yaml"

wildcard_constraints:
    datasource = "Geuvadis|GTEx",
    population = "CEU|FIN|GBR|TSI|YRI|EUR",
    geuvadis_sid = "[A-Z]{2}[0-9]{5,6}",
    chrom = "chr[0-9]{1,2}",
    window = "[0-9]+"

rule RunLeafcutter:
    input:
        junc_files = 'resources/{datasource}/{population}/juncs',
    output:
        counts = 'results/pheno/leafcutter/{datasource}/{population}/leafcutter_perind.counts.gz'
    params:
        run_dir    = 'results/pheno/leafcutter/{datasource}/{population}',
        out_prefix = 'leafcutter', # file name prefix, not directory
        py_script = 'workflow/submodules/leafcutter/clustering/leafcutter_cluster_regtools.py'
    conda: 'leafcutter'
    threads: 1
    resources: cpu=1, time=2100, mem_mb=25000
    shell:
        '''
        python {params.py_script} \
            -j <(realpath {input.junc_files}/*autosomes.junc) \
            -o {params.out_prefix} \
            -r {params.run_dir} 
        
        '''




rule FixSampleNames:
    input: 
        counts = 'results/pheno/leafcutter/{datasource}/{population}/leafcutter_perind.counts.gz'
    output:
        counts = temp('results/pheno/leafcutter/{datasource}/{population}/leafcutter_perind_fixednames.counts.gz')
    threads: 1
    run:
        import gzip

        fout = output.counts
        with gzip.open(input.counts) as fin:
            with gzip.open(output.counts, 'wt') as fout:
                i = 0
                for line in fin:
                    line = line.decode('utf-8').strip().split(' ')
                    line = [x.strip() for x in line]

                    if i == 0:
                        
                        # from second elements, extract sample ids delimited by '.'
                        samples = [x.split('.')[0] for x in line[1:]]

                        # cat message with sample ids
                        print(f'Fixing column names of {len(samples)} samples:\n\n{", ".join(samples)}')

                        # rename header using first element of line and samples
                        line = [line[0]] + samples

                        # join back to string with ' ' and write out
                        fout.write(' '.join(line) + '\n')

                    else:
                        fout.write(' '.join(line) + '\n')

                    i += 1


def getExtractCountsInput(wildcards):
    if wildcards.datasource == 'GTEx':
        counts = os.path.join('results/pheno/leafcutter', wildcards.datasource, 
                wildcards.population, 'leafcutter_perind_fixednames.counts.gz')
        vcf = config['annotation']['GTEx']['VCF_v7_HG38']
    elif wildcards.datasource == 'Geuvadis':
        counts = os.path.join('results/pheno/leafcutter', wildcards.datasource, 
                wildcards.population, 'leafcutter_perind_fixednames.counts.gz')
        vcf = '/project2/yangili1/zpmu/1kg_b38/CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz'
        # only need the vcf header to get all available sample names, thus choose a smaller vcf
    return {'counts': counts, 'vcf': vcf}

rule ExtractCounts:
    input: unpack(getExtractCountsInput)
    output:
        counts = 'results/pheno/leafcutter/{datasource}/{population}/leafcutter_perind_GenoMatched.counts.gz',
        indivs = 'results/pheno/leafcutter/{datasource}/{population}/individuals.txt'
    run:
        import subprocess
        import gzip
        vcf_header = f"bcftools view -h {input.vcf} | awk '$1 ~ /#CHROM/' "
        vcf_header = subprocess.run(vcf_header, shell=True, capture_output=True).stdout.decode('utf-8').strip().split()
        geno_inds = vcf_header[9:] # available genotype individual IDs

        fout_counts = gzip.open(output.counts, 'wt')
        with gzip.open(input.counts) as f:
            i = 0
            for ln in f:
                lnsplit = ln.decode().strip().split()
                if i == 0: # header
                    keep1 = [lnsplit.index(l) for l in lnsplit if l in ['chrom'] + geno_inds ] # cols to keep
                    keep2 = [lnsplit.index(l) for l in lnsplit if l in geno_inds ] # cols to keep, no 'chrom'

                    # write individual ids to a text file
                    with open(output.indivs, 'w') as sf: # write individual ids
                        indvs = [lnsplit[ix] for ix in keep2]
                        print(f'Removing these samples not in genotype vcf:\n{[x for x in lnsplit[1:] if x not in indvs]}')
                        sf.writelines([x + '\n' for x in indvs])
                    
                    # write header line for counts
                    buf = ' '.join([lnsplit[ix] for ix in keep1]) + '\n' # header line
                    fout_counts.write(buf)
                    print(f'Finished writing sample ids into: {output.indivs}')
                else:
                    buf = ' '.join([lnsplit[ix] for ix in keep1]) + '\n' # header line
                    fout_counts.write(buf)
                
                i += 1
                
        fout_counts.close()
        print(f'Finished writing counts file into: {output.counts}') 




rule PrepPhenoBed:
    message: '### Prepare phenotype bed file for qtltools (rank normalized)'
    input:
        counts = 'results/pheno/leafcutter/{datasource}/{population}/leafcutter_perind_GenoMatched.counts.gz',
        anno = config['annotation']['gencode']
    output: 'results/pheno/leafcutter/{datasource}/{population}/pheno.chr1.bed.gz'
    params:
        Rscript = 'workflow/scripts/prepPhenoBed.R',
        outprefix = 'results/pheno/leafcutter/{datasource}/{population}/pheno',
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
        
        ls {output}
        '''


rule PhenotypePCA:
    message: '### Run PCA on phenotype with permutation'
    input: rules.PrepPhenoBed.output
    output: 'results/pheno/leafcutter/{datasource}/{population}/{chrom}.pca'
    params:
        rscript = 'workflow/scripts/PermuteAndPCA.R',
        inputfile = 'results/pheno/leafcutter/{datasource}/{population}/pheno.{chrom}.bed.gz'
    shell:
        '''
        ls {input}
        Rscript {params.rscript} {params.inputfile} {output} 
        '''



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
    message: '### Extract genotype vcf for phenotype samples - {wildcards.datasource}:{wildcards.population}:{wildcards.chrom}'
    input: 
        vcf = getExtractGenotypeInput,
        sample_file = rules.ExtractCounts.output.indivs
    output: 'results/geno/{datasource}/{population}/{chrom}.vcf.gz'
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
    input:  rules.ExtractGenotypeVCF.output
    output: 'results/geno/{datasource}/{population}/{chrom}.pca'
    params:
        out_prefix = 'results/geno/{datasource}/{population}/{chrom}'
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
    output: 'results/pheno/leafcutter/{datasource}/{population}/{chrom}_CovMatrix.txt'
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
        vcf = 'results/geno/{datasource}/{population}/{chrom}.vcf.gz',
        bed = 'results/pheno/leafcutter/{datasource}/{population}/pheno.{chrom}.bed.gz',
        cov = 'results/pheno/leafcutter/{datasource}/{population}/{chrom}_CovMatrix.txt'
    output: temp('results/qtl/leafcutter/{datasource}/{population}/cis_{window}/perm/{chrom}.txt')
    params:
        cis_window = '{window}'
    resources: cpu = 1, mem = 12000, time = 1000
    shell:
        '''
        QTLtools cis \
            --seed 123 --silent \
            --vcf {input.vcf} --bed {input.bed} --cov {input.cov}  --out {output} \
            --window {params.cis_window} \
            --permute 1000 --region {wildcards.chrom}
        '''


rule AddQvalueToPermutationPass:
    message: '### Add qvalue to the output of QTLtools permutation pass'
    input:  'results/qtl/leafcutter/{datasource}/{population}/cis_{window}/perm/{chrom}.txt'
    output: 'results/qtl/leafcutter/{datasource}/{population}/cis_{window}/perm/{chrom}.addQval.txt.gz'
    params:
        rscript = 'workflow/scripts/AddQvalueToQTLtoolsOutput.R'
    shell:
        '''
            txt=$(echo {output} | sed -E 's/.gz//')
            Rscript {params.rscript} {input} $txt
            bgzip $txt
        '''


