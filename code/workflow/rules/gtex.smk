'''
NOTE:
    - `SMTS` and the `SMTSD` fields in the metadata has white spaces. White spaces
      removed in dataframe. When using tissue wildcards, just remove white spaces.


See details about GTEx tissues and tissue subcategories here:
https://daiuc.github.io/SpliFi/analysis/2023-11-17-GTEx-tissues-summary.pub.html



'''



rule MungeGTExJuncs:
    message: '# Munge GTEx junc file'
    input: 'resources/GTEx/juncs/GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct.gz'
    output: 
        'resources/GTEx/juncs/munged_juncs.txt.gz'
    params: 
        anno = config['annotation']['gencode_v26_genes'], # used to infer strand using gene_id; input file doesn't have strand info
        py_script = 'workflow/scripts/munge_gtex_juncs.py'
    threads: 1
    resources: cpu = 1, mem_mb = 15000, time = 2100
    shell: 
        '''
        python {params.py_script} -I {input} -O {output} -G {params.anno}
        '''

rule SplitGTExSampleIDList:
    '''
    This step is necessary purely for managing workflow. There 
    are 17382 samples, each sample times about 100s to split. 
    This rule split these samples into 70 groups, saving 
    sample IDs in 70 files. The next rule the parallize using
    70 compute nodes to split all samples within each file.
    '''
    input: 'resources/GTEx/juncs/SampleIDs.txt'
    output: touch('resources/GTEx/juncs/SampleIDs.split.done')
    params:
        Num_Groups = 70, # number of parts to split into
        out_prefix = 'resources/GTEx/juncs/SampleIDs.part_'
    run:
        from math import ceil
        with open(input[0]) as f:
            SAMIDS = [c.strip() for c in f.readlines()]
            R = params.Num_Groups # number of split parts
            L = ceil(len(SAMIDS)/R) # number of samples per part
            for i in range(R):
                if len(SAMIDS) >= i*L + L:
                    WriteIDs = SAMIDS[i*L : i*L + L]
                else:
                    WriteIDs = SAMIDS[i*L : len(SAMIDS)]
                WriteFN = params.out_prefix + str(i) + '.txt'
                with open(WriteFN, 'w') as fo:
                    for s in WriteIDs:
                        fo.write(s + '\n')
        

rule SplitGTExJuncs:
    message: '# Split munged GTEx junc files by sample'
    input: 
        counts = 'resources/GTEx/juncs/munged_juncs.txt.gz',
        sample_file = 'resources/GTEx/juncs/SampleIDs.part_{GTEx_junc_part}.txt'
    output: touch('resources/GTEx/juncs/split_juncs/part_{GTEx_junc_part}.done')
    wildcard_constraints:
        GTEx_junc_part = '[0-9]{1,2}'
    params:
        py_script = 'workflow/scripts/split_gtex_juncs.py',
        out_prefix = 'resources/GTEx/juncs/split_juncs/' 
    threads: 1
    resources: cpu = 1, mem_mb = 15000, time = 2000
    shell: 
        '''
        # for zsh IFS=$'\\n' arr=($(cat path/to/file))
        # for bash
        readarray -t SAMPLE_IDS < {input.sample_file}
        for s in ${{SAMPLE_IDS[@]}}; do
            OUT_FILE={params.out_prefix}${{s}}.tsv.gz
            echo running: python {params.py_script} -I {input.counts} -O $OUT_FILE -S $s
            python {params.py_script} -I {input.counts} -O $OUT_FILE -S $s
        done
        '''




# ----------------------------------------------------------------------------------------




# aggregate junction file by subject id
# when running with sub-tissue types, there is no combining, since each tissue type
# only has 1 tissue per individual (subject id)
rule AggregateJuncBySubjGtex:
    '''
    Aggregate junction files by subject id. Because in the GTEx dataset,
    each subject ID corresponds to a person, but the original junctil files
    are organized per sample id. Each Subject ID can have 1 or more sample IDs.
    So this step creates a single junction file for each individual.
    '''
    output: 
        out_flag = temp(touch('resources/GTEx/juncs/groupped_juncs/{tissue}/uncompressed/groupjunc.done')),
        out_dir = temp(directory('resources/GTEx/juncs/groupped_juncs/{tissue}/uncompressed'))
    params:
        inputDir = 'resources/GTEx/juncs/split_juncs',
        metadata = config['Dataset']['GTEx']['Junc_meta'],
        pyscript = 'workflow/scripts/aggJuncs.py'
    log: 'logs/AggregateJunBySubjGtex/{tissue}.log'
    shell:
        '''
        python {params.pyscript} --metadata {params.metadata} \
            --inDir {params.inputDir} --tissue "{wildcards.tissue}" \
            --outDir {output.out_dir} &> {log}
        
        '''
                

rule CompressGrouppedJuncs:
    input: rules.AggregateJuncBySubjGtex.output.out_dir
    output: 
        flag = temp(touch('resources/GTEx/juncs/groupped_juncs/{tissue}/compressed/done')),
        out_dir = temp(directory('resources/GTEx/juncs/groupped_juncs/{tissue}/compressed'))
    resources: cpu=8, time=1200, mem_mb=25000
    threads: 8
    log: 'logs/CompressGrouppedJuncs/{tissue}.log'
    shell:
        '''
        # module load parallel

        # Find all .tsv files in the input folder and use parallel to gzip them
        find "{input}" -type f -name "*.tsv" -print0 | \
            parallel -0 -j {threads} "gzip -c {{}} > {output.out_dir}/{{/.}}.tsv.gz"
        
        ls 

        '''


rule ConvertCoordinates:
    '''
    I found out that the junction coordinates in the GTEx junction files are
    VCF-like coordinates, which is not compatible with the coordinates in the
    annotation file. This rule converts the coordinates in the junction files
    to the annotation file coordinates - aka BED like coordinates.
    '''
    input: 
        flag = rules.CompressGrouppedJuncs.output.flag,
        in_Dir = rules.CompressGrouppedJuncs.output.out_dir
    output: 
        flag = touch('resources/GTEx/juncs/groupped_juncs/{tissue}/converted/done'),
        out_dir = directory('resources/GTEx/juncs/groupped_juncs/{tissue}/converted')
    params:
        input_dir = 'resources/GTEx/juncs/groupped_juncs/{tissue}/compressed'
    resources: cpu=8, time=1200, mem_mb=25000
    threads: 8
    log: 'logs/ConvertCoordinates/{tissue}.log'
    shell:
        '''
        # Use awk to subtract the second column by 1
        # module load parallel

        echo converting junction coordinates to BED-like coordinates ... > {log}
        
        subtractStart() {{
            in_f=$1
            ou_f=$2
            zcat $in_f | \
                awk -v OFS='\t' '{{print $1, $2-1, $3, $4, $5, $6}}' | \
                gzip -c > $ou_f
            
            echo "Done converting $in_f -> $ou_f" >> {log}
        }}

        export -f subtractStart

        find {input.in_Dir} -type f -name "*.gz" -print0 | \
            parallel -0 -j {threads} "subtractStart {{}} {output.out_dir}/{{/.}}.gz"

        '''
    
rule RenameJuncFilesGtex:
    message: '### Rename GTEx junction files to include tissue name and SUBJID'
    input: rules.ConvertCoordinates.output.flag
    output: touch('resources/GTEx/juncs/all49tissues/{tissue}.done')
    params:
        source_dir = 'resources/GTEx/juncs/groupped_juncs/{tissue}/converted',
        dest_prefix = 'resources/GTEx/juncs/all49tissues/{tissue}'
    log: 'logs/RenameJuncFilesGtex/{tissue}.log'
    shell:
        '''
        source_files=({params.source_dir}/*.tsv.gz)
        for f in ${{source_files[@]}}; do
            subjID=$(basename -s .tsv.gz $f)
            dest_file={params.dest_prefix}.${{subjID}}.tsv.gz

            cp $f $dest_file && zcat $dest_file | wc -l &>> {log} &

            if (( $(jobs | wc -l) > 20 )); then
                echo "wait for previous jobs to finish" &>> {log}
                wait
            fi
            
        done
        echo "All done $(date)" &>> {log}
        '''



# ----------------------------------------------------------------------------------------



# optional rule used to speed up the workflow, intron clusters may run once
# and reuse as needed.
rule ClusterIntronsGtex:
    message: '### Make intron clusters using any GTEx samples'
    input:
        junc_files_flag = 'resources/GTEx/juncs/groupped_juncs/{tissue}/converted/done',
    output: 
        pooled   = 'results/pheno/noisy/GTEx/{tissue}/leafcutter_pooled',
        clusters = 'results/pheno/noisy/GTEx/{tissue}/leafcutter_refined_noisy',
        lowusage = 'results/pheno/noisy/GTEx/{tissue}/leafcutter_lowusage_introns', # intermediate
        refined  = 'results/pheno/noisy/GTEx/{tissue}/leafcutter_refined' # intermediate
    params:
        run_dir    = 'results/pheno/noisy/GTEx/{tissue}',
        out_prefix = 'leafcutter',
        junc_files = 'resources/GTEx/juncs/groupped_juncs/{tissue}/converted',
        py_script  = 'workflow/submodules/leafcutter2/scripts/leafcutter_make_clusters.py'
    log: 'logs/ClusterIntronsGtex/{tissue}.log'
    shell:
        '''
        python {params.py_script} \
            -r {params.run_dir} \
            -o {params.out_prefix} \
            -j <(realpath {params.junc_files}/*.tsv.gz) &> {log}
        
        ls -lah {output.pooled} {output.clusters} {output.lowusage} {output.refined} &>> {log}
        
        '''


rule Leafcutter2Gtex:
    message:'### Run leafcutter2 on GTEx samples'
    input: 
        junc_files_flag = 'resources/GTEx/juncs/groupped_juncs/{tissue}/converted/done',
        junc_files = 'resources/GTEx/juncs/groupped_juncs/{tissue}/converted',
        intron_clusters = 'results/pheno/noisy/GTEx/{tissue}/leafcutter_refined_noisy',
    output:
        perind_noise_by_intron = 'results/pheno/noisy/GTEx/{tissue}/leafcutter_perind.counts.noise_by_intron.gz'
    params:
        py_script  = 'workflow/submodules/leafcutter2/scripts/leafcutter2_regtools.py',
        run_dir    = 'results/pheno/noisy/GTEx/{tissue}',
        out_prefix = 'leafcutter', # note do not include parent dir
        gtf = config['annotation']['gtf']['v43'],
        genome = config['genome38'],
        pre_clustered = '-c results/pheno/noisy/GTEx/{tissue}/leafcutter_refined_noisy', 
        other_params = '-k' # not keeping constitutive introns
    threads: 1
    resources: cpu=1, time=2100, mem_mb=25000
    log: 'logs/Leafcutter2Gtex/{tissue}.log'
    shell:
        '''
        python {params.py_script} \
            -j <(realpath {input.junc_files}/*.tsv.gz) \
            -r {params.run_dir} \
            -o {params.out_prefix} \
            -A {params.gtf} \
            -G {params.genome} \
            {params.pre_clustered} {params.other_params} &> {log}

        ls {output.perind_noise_by_intron} &>> {log}
        

        '''


use rule Leafcutter2Gtex as Leafcutter2Gtex_wConst with:
    message: '### run leafcutter2 on GTEx samples with constitutive introns'
    input: 
        junc_files_flag = 'resources/GTEx/juncs/groupped_juncs/{tissue}/converted/done',
        junc_files = 'resources/GTEx/juncs/groupped_juncs/{tissue}/converted',
    output:
        perind_noise_by_intron = 'results/pheno/noisy/GTEx/{tissue}/wConst_perind.constcounts.noise_by_intron.gz'
    params:
        py_script  = 'workflow/submodules/leafcutter2/scripts/leafcutter2_regtools.py',
        run_dir    = 'results/pheno/noisy/GTEx/{tissue}',
        pre_clustered = '-c results/ds/GTEx/all49tissues_refined_noisy',
        out_prefix = 'wConst',
        gtf = config['annotation']['gtf']['v43'],
        genome = config['genome38'],
        other_params = '-k --includeconst' # include constitutive introns
    log: 'logs/Leafcutter2Gtex_wConst/{tissue}.log'





# ----------------------------------------------------------------------------------------
#           Cluster introns using all GTEx tissues at once
#         Then use the clusters to run leafcutter2 on specific tissues
#         This enables differential splicing analysis across tissues
# ----------------------------------------------------------------------------------------


# optional rule used to speed up the workflow, intron clusters may run once
# and reuse as needed.
rule ClusterIntronsGtexAllTissues:
    message: '### Make intron clusters using all (49) GTEx tissues'
    output: 
        pooled   = 'results/ds/GTEx/all49tissues_pooled',
        clusters = 'results/ds/GTEx/all49tissues_refined_noisy',
        lowusage = 'results/ds/GTEx/all49tissues_lowusage_introns', # intermediate
        refined  = 'results/ds/GTEx/all49tissues_refined' # intermediate
    params:
        run_dir    = 'results/ds/GTEx',
        out_prefix = 'all49tissues',
        junc_files = 'resources/GTEx/juncs/all49tissues',
        py_script  = 'workflow/submodules/leafcutter2/scripts/leafcutter_make_clusters.py'
    log: 'logs/ClusterIntronsGtexAllTissues/all49tissues.log'
    resources: cpu=1, time=2100, mem_mb=35000
    shell:
        '''
        python {params.py_script} \
            -r {params.run_dir} \
            -o {params.out_prefix} \
            -j <(ls {params.junc_files}/*.tsv.gz) &> {log}
        
        ls -lah {output.pooled} {output.clusters} {output.lowusage} {output.refined} &>> {log}
        
        '''


rule LeafcutterForDSGtex:
    message: '### Run leafcutter2 on GTEx samples for differential splicing analysis'
    input:
        pre_clusters = 'results/ds/GTEx/all49tissues_refined_noisy',
        tissue1_flag = 'resources/GTEx/juncs/all49tissues/{ds_tissue_1}.done',
        tissue2_flag = 'resources/GTEx/juncs/all49tissues/{ds_tissue_2}.done',
    output:
        ds_numers = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds_perind_numers.counts.noise_by_intron.gz',
        # ds_counts_lf1 is necessary for leafcutter_ds.R script
        ds_numers_lf1 = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds_perind_numers.counts.noise_by_intron.lf1.gz',
        ds_sample_group = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds_sample_group.txt'
    params:
        run_dir = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}',
        out_prefix = 'ds',
        tissue1_juncs = 'resources/GTEx/juncs/all49tissues/{ds_tissue_1}',
        tissue2_juncs = 'resources/GTEx/juncs/all49tissues/{ds_tissue_2}',
        NSamples = 50, # select this number of samples per tissue type
        pre_clustered = '-c results/ds/GTEx/all49tissues_refined_noisy',
        gtf = config['annotation']['gtf']['v43'],
        genome = config['genome38'],
        other_params = '-k ', # not keeping constitutive introns
        py_script  = 'workflow/submodules/leafcutter2/scripts/leafcutter2_regtools.py',
        py_script2 = 'workflow/scripts/makeSampleGroupFileForDifferentialSplicing.py'
    log: 'logs/LeafcutterForDSGtex/{ds_tissue_2}_v_{ds_tissue_1}.log'
    resources: cpu = 1, time = 2100, mem_mb = 25000
    shell:
        '''
        # run leafcutter2 for differential analysis
        python {params.py_script} \
            -j <(cat <(ls {params.tissue1_juncs}*tsv.gz | head -{params.NSamples}) <(ls {params.tissue2_juncs}*tsv.gz | head -{params.NSamples})) \
            -r {params.run_dir} \
            -o {params.out_prefix} \
            -A {params.gtf} \
            -G {params.genome} \
            {params.pre_clustered} {params.other_params} &> {log}

        # make sample group file for differential splicing analysis
        python {params.py_script2} -i {output.ds_numers} -o {output.ds_numers_lf1} -s {output.ds_sample_group} &>> {log}

        '''


rule RunLeafcutterDiffSplicingGtex:
    input:
        ds_numers_lf1 = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds_perind_numers.counts.noise_by_intron.lf1.gz',
        ds_sample_group = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds_sample_group.txt'
    output:
        flag = touch('results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/done')
        # produces two files:
        # 1. {outprefix}_effect_sizes.txt
        # 2. {outprefix}_manual_ds_cluster_significance.txt
    params:
        Rscript = 'workflow/submodules/leafcutter/scripts/leafcutter_ds.R', 
        outprefix = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds', # note you need to include path!
        MIN_SAMPLES_PER_INTRON = 5,
        MIN_SAMPLES_PER_GROUP = 3,
        MIN_COVERAGE = 5
    resources: cpu = 4, mem_mb = 30000, time = 2100
    threads: 4
    conda: 'leafcutter1'
    log: 'logs/RunLeafcutterDiffSplicingGtex/{ds_tissue_2}_v_{ds_tissue_1}.log'
    shell:
        '''
        {params.Rscript} --num_threads {threads} \
            --output_prefix {params.outprefix} \
            --min_samples_per_intron={params.MIN_SAMPLES_PER_INTRON} \
            --min_samples_per_group={params.MIN_SAMPLES_PER_GROUP} \
            --min_coverage={params.MIN_COVERAGE} \
            {input.ds_numers_lf1} {input.ds_sample_group} &> {log}


        '''

# NOTE to make leafcutter's leafcutter_ds.R script work, 
# must modify the _perind_numers.counts.noise_by_intron.gz file
# First row should not have the 'chrom' first column
# first columns can only be like `chr1:827775:829002:clu_1_+` because 
# the R function only expect to split by ":" into 4 columns.



# ----------------------------------------------------------------------------------------
# ad hoc

rule adhoc_test_ds_step1: # testing differential splicing with neg control
    input: 
        ds_numers_lf1 = 'results/ds/GTEx/Brain-Cerebellum_v_Liver/ds_perind_numers.counts.noise_by_intron.lf1.gz',
    output:
        neg_control_numers = 'results/ds/GTEx/ds_test/BC_v_Liver/perind.numers.gz',
        neg_control_groups = 'results/ds/GTEx/ds_test/BC_v_Liver/sample_group.txt'
    run:
        import gzip
        outf1 = gzip.open(output.neg_control_numers, 'wt')
        outf2 = open(output.neg_control_groups, 'w')
        with gzip.open(input.ds_numers_lf1, 'rt') as f:
            i = 0
            for ln in f:
                if i == 0:
                    header = ln.split()
                    cols = ([x.replace('Brain-Cerebellum', 'BC-group1') for x in header[0:100]] +
                            [x.replace('Brain-Cerebellum', 'BC-group2') for x in header[100:200]])
                    groups = ['BC-group1' for x in header[0:100]] + ['BC-group2' for x in header[100:200]]
                    for c,g in zip(cols, groups):
                        outf2.write(f'{c} {g}\n')
                    outf1.write(' '.join(cols) + '\n')
                if i > 0:
                    outln = ln.split()[:201]
                    outf1.write(' '.join(outln) + '\n')
                i += 1

        outf1.close()
        outf2.close()


use rule RunLeafcutterDiffSplicingGtex as adhoc_test_ds_step2 with:
    input:
        ds_numers_lf1 = 'results/ds/GTEx/ds_test/BC_v_Liver/perind.numers.gz',
        ds_sample_group = 'results/ds/GTEx/ds_test/BC_v_Liver/sample_group.txt'
    output:
        flag = touch('results/ds/GTEx/ds_test/BC_v_Liver/ds.done')
        # produces two files:
        # 1. {outprefix}_effect_sizes.txt
        # 2. {outprefix}_manual_ds_cluster_significance.txt
    params:
        Rscript = 'workflow/submodules/leafcutter/scripts/leafcutter_ds.R', 
        outprefix = 'results/ds/GTEx/ds_test/BC_v_Liver/ds', # note you need to include path!
        MIN_SAMPLES_PER_INTRON = 5,
        MIN_SAMPLES_PER_GROUP = 3,
        MIN_COVERAGE = 5
    log: 'results/ds/GTEx/ds_test/BC_v_Liver/log'



## -----------------------------------------------------------------------------
##   GTEx expression data
## -----------------------------------------------------------------------------

rule ExtractGTExGeneExpression:
    input: 
      tpm = 'resources/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz',
      cnt = 'resources/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz'
    output: 
      tpm = 'resources/GTEx/expression/{tissue}_gene_tpm.tsv.gz',
      cnt = 'resources/GTEx/expression/{tissue}_gene_reads.tsv.gz',
    params: 
        py_script = 'workflow/scripts/extract_gtex_gene_expression.py',
        junc_meta = config['Dataset']['GTEx']['Junc_meta']
    log: 'logs/ExtractGTExGeneExpression/{tissue}.log'
    shell:
        '''
        echo extracting tpm from {input.tpm} ...&> {log}
        python {params.py_script} \
            -I {input.tpm} \
            -M {params.junc_meta} \
            -O {output.tpm} \
            -T {wildcards.tissue} &>> {log}

        echo extracting raw counts from {input.cnt} ...
        python {params.py_script} \
            -I {input.cnt} \
            -M {params.junc_meta} \
            -O {output.cnt} \
            -T {wildcards.tissue} &>> {log}

        echo "Number of lines in {output.tpm} : $(zcat {output.tpm} | wc -l)" &>> {log}
        echo "Number of lines in {output.cnt} : $(zcat {output.cnt} | wc -l)" &>> {log}
        '''

# prepare GTEX dge data
rule PrepareGTExDGE:
    input: 
        cnt1 = 'resources/GTEx/expression/{dge_tissue1}_gene_reads.tsv.gz',
        cnt2 = 'resources/GTEx/expression/{dge_tissue2}_gene_reads.tsv.gz',
    output: 
        #NOTE: tissue2 is intended to be numerator, and tissue 1 denominator, in subsequent dge step
        cnt = 'results/dge/GTEx/{dge_tissue2}_v_{dge_tissue1}_counts.tsv',
        coldata = 'results/dge/GTEx/{dge_tissue2}_v_{dge_tissue1}_coldata.tsv',
    params:
        R_script = 'workflow/scripts/prepare_GTEx_dge.R',
        n_samples = 50, # number of samples to use for each tissue
        outdir = 'results/dge/GTEx'
    log: 'logs/PrepareGTExDGE/{dge_tissue2}_v_{dge_tissue1}.log'
    shell:
        '''
        Rscript {params.R_script} {input.cnt1} {input.cnt2} {params.n_samples} {params.outdir} &> {log}
        ls {output.cnt} {output.coldata} &>> {log}
        '''

# run dge
rule DgeGtex:
  input: 
    cnt = 'results/dge/GTEx/{dge_tissue2}_v_{dge_tissue1}_counts.tsv',
    coldata = 'results/dge/GTEx/{dge_tissue2}_v_{dge_tissue1}_coldata.tsv'
  output:
    dge = 'results/dge/GTEx/{dge_tissue2}_v_{dge_tissue1}_dge_genes.tsv',
  params:
    R_script = 'workflow/scripts/dge.R',
    outprefix = 'results/dge/GTEx/{dge_tissue2}_v_{dge_tissue1}',
    min_reads = 10,
    min_samples = 10,
    min_fdr = 0.1,
    min_log2fc = 1,
  log: 'logs/DgeGtex/{dge_tissue2}_v_{dge_tissue1}.log'
  shell:
    '''
    Rscript {params.R_script} \
            {input.cnt} {input.coldata} {params.outprefix} \
            {params.min_reads} {params.min_samples} {params.min_fdr} {params.min_log2fc} &> {log}
    ls {output.dge} &>> {log}

    '''







## -----------------------------------------------------------------------------
##   plot sashimi
## -----------------------------------------------------------------------------

def get_bedgraph_input(wildcards):
    import glob

    bedgraph_dir = '/project2/yangili1/GTEx_v8/bedGraph'
    tissueTrans = GTEX_BEDGRAPH_TISSUES.get(wildcards.tissue)
    bedFiles = glob.glob(f'{bedgraph_dir}/{tissueTrans}/*.bed.gz') # abs paths

    return bedFiles

rule BedgraphToBW:
    '''
    Note the wc for tissue here is the original tissue name, not transformed
    '''
    input: get_bedgraph_input
    output: touch('resources/GTEx/BigWig/{tissue}/done')
    params:
        outPrefix = 'resources/GTEx/BigWig/{tissue}',
        chromSizes = 'resources/hg38_w_chrEBV.chrom.sizes'
    threads: 8
    resources: cpu=8, mem_mb=25000, time=1200
    shell:
        '''
        bedFiles="{input}"

        numJobs=0
        maxJobs={threads}

        for f in ${{bedFiles[@]}}; do
            tmpf={params.outPrefix}/$(basename "$f" .bed.gz).tmp.bed
            outf={params.outPrefix}/$(basename $f .bed.gz).bw

            if [[ $numJobs -le $maxJobs ]]; then
                echo bedGraphToBigWig $f ...
                bgzip -f -d -c $f > $tmpf && bedGraphToBigWig $tmpf {params.chromSizes} $outf && rm $tmpf &
                numJobs=$((numJobs + 1))
            else
                wait
                numJobs=$((numJobs - 1))
                
                echo bedGraphToBigWig $f ...
                bgzip -f -d -c $f > $tmpf && bedGraphToBigWig $tmpf {params.chromSizes} $outf && rm $tmpf &
                numJobs=$((numJobs + 1))
            fi
        done

        wait
        echo "All done!"

        '''


