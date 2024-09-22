'''
Differential splicing analysis &
Differential expression analysis &
related analysis & data preps & plots

'''

N_DIFFER = 80 # number of samples to choose for differnetial splicing and expression analysis


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
        NSamples = N_DIFFER, # select this number of samples per tissue type
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
    message: """Run differential splicing analysis using leafcutter1's leafcutter_ds.R script"""
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
# NOTE: need to rerun this step and below because right now samples
# are not matched with ds samples! 
# Should use the ds_sample_group.txt file to get matching samples for each tissue
rule PrepareGTExDGE:
    input: 
        cnt1 = 'resources/GTEx/expression/{dge_tissue1}_gene_reads.tsv.gz',
        cnt2 = 'resources/GTEx/expression/{dge_tissue2}_gene_reads.tsv.gz',
        ds_samples = 'results/ds/GTEx/{dge_tissue2}_v_{dge_tissue1}/ds_sample_group.txt'
    output: 
        #NOTE: tissue2 is intended to be numerator, and tissue 1 denominator, in subsequent dge step
        cnt = 'results/dge/GTEx/{dge_tissue2}_v_{dge_tissue1}_counts.tsv',
        coldata = 'results/dge/GTEx/{dge_tissue2}_v_{dge_tissue1}_coldata.tsv',
    params:
        R_script = 'workflow/scripts/prepare_GTEx_dge.R',
        outdir = 'results/dge/GTEx'
    log: 'logs/PrepareGTExDGE/{dge_tissue2}_v_{dge_tissue1}.log'
    shell:
        '''
        Rscript {params.R_script} {input.cnt1} {input.cnt2} {input.ds_samples} {params.outdir} &> {log}
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
  log: 'logs/DgeGtex/{dge_tissue2}_v_{dge_tissue1}.log'
  shell:
    '''
    Rscript {params.R_script} \
            {input.cnt} {input.coldata} {params.outprefix} \
            {params.min_reads} {params.min_samples}  &> {log}
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


#-----------------------------------------------------------------------------------------
#   plot sashimi for differential splicing
#-----------------------------------------------------------------------------------------

rule getIntronsForSashimi:
    input:
        rds = '../data/ds_v_dge/{ds_tissue_2}_v_{ds_tissue_1}_data.rds' # separate smk rules to make rds files 
    output: 'plots/sashimi/ds/{ds_tissue_2}_v_{ds_tissue_1}/introns.txt'
    params:
        py_script = 'workflow/scripts/getIntrons.R',
        contrast = '{ds_tissue_2}_v_{ds_tissue_1}',
        minDeltaPsi = 0.2,
        FDR = 0.001,
        minL2FC = 0.58,
    log: 'logs/getIntronsForSashimi/{ds_tissue_2}_v_{ds_tissue_1}.log'
    shell:
        '''
        Rscript {params.py_script} {input.rds} {params.contrast} {output} {params.minDeltaPsi} {params.FDR} {params.minL2FC} &> {log}
        '''


def getPlotIntrons(wildcards):
    with open(f'plots/sashimi/ds/{wildcards.ds_tissue_2}_v_{wildcards.ds_tissue_1}/introns.txt') as f:
        introns = [x.strip() for x in f.readlines()]
    return introns

rule PrepSashimiDsGtex:
    message: '### Prepare sashimi plots for differential splicing'
    input: 
        ds_sample_group = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds_sample_group.txt',
        ds_effect_sizes = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds_effect_sizes.txt',
        ds_intron_count = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds_perind_numers.counts.noise_by_intron.gz',
        introns = 'plots/sashimi/ds/{ds_tissue_2}_v_{ds_tissue_1}/introns.txt',
    output:
        flag = touch('plots/sashimi/ds/{ds_tissue_2}_v_{ds_tissue_1}/prep.done')
        # flag = touch('plots/sashimi/ds/{ds_tissue_2}_v_{ds_tissue_1}/{plotIntron}.prep.done')
    params:
        contrast = '{ds_tissue_2}_v_{ds_tissue_1}',
        outDir = 'plots/sashimi/ds/{ds_tissue_2}_v_{ds_tissue_1}',
        bwPrefix = 'resources/GTEx/BigWig',
        plotIntrons = getPlotIntrons,
        iniTemplate = 'config/template-sashimi-diffsplice.ini',
        shellTemplate = 'config/template-plot-sashimi-cmd.sh',
        pyscript = 'workflow/scripts/prepSashimi.py',
    threads: 4
    resources: cpu=1, mem_mb=25000, time=1200
    group: 'PrepSashimiDsGtex'
    log: 'logs/PrepSashimiDsGtex/{ds_tissue_2}_v_{ds_tissue_1}.log'
    shell:
        '''
        #module load parallel

        introns="{params.plotIntrons}"

        cmd="python {params.pyscript}" 
        cmd+=" --contrast {params.contrast} "
        cmd+=" --outDir {params.outDir} "
        cmd+=" --bwPrefix {params.bwPrefix} "
        cmd+=" --dsSampleGroupFile {input.ds_sample_group} "
        cmd+=" --dsEffectFile {input.ds_effect_sizes} "
        cmd+=" --intronCountsFile {input.ds_intron_count} "
        cmd+=" --plotIntron {{}} "
        cmd+=" --iniTemplate {params.iniTemplate} "
        cmd+=" --plotShellTemplate {params.shellTemplate} "
                                 
        parallel -j {threads} "$cmd" ::: $introns &> {log}

        '''

def getPLotSashimiDsGtexParams(wildcards):
    w = wildcards
    folder = f'plots/sashimi/ds/{w.ds_tissue_2}_v_{w.ds_tissue_1}'
    clu = w.plotIntron.split(':')[-1]

    try:
        shellFiles = [os.path.basename(x) for x in glob.glob(f'{folder}/*.sh')]
        shell = [x for x in shellFiles if clu in x][0]
    except:
        shell = None

    return f'{shell}'

rule PlotSashimiDsGtex:
    input: 'plots/sashimi/ds/{ds_tissue_2}_v_{ds_tissue_1}/prep.done'
    output: touch('plots/sashimi/ds/{ds_tissue_2}_v_{ds_tissue_1}/plot.done')
    params:
        folder = 'plots/sashimi/ds/{ds_tissue_2}_v_{ds_tissue_1}',
        # shell = getPLotSashimiDsGtexParams
    conda: 'pygenometracks'
    group: 'PrepSashimiDsGtex'
    threads: 4
    log: 'logs/PlotSashimiDsGtex/{ds_tissue_2}_v_{ds_tissue_1}.log'
    shell: 
        '''
        # module load parallel

        log=$(realpath {log})
        echo "plot sashimi plots" > $log
        cd {params.folder}
        shellFiles=$(ls *.sh)
        parallel -j {threads} "sh {{}}" ::: $shellFiles &>> $log

        #run in login node after: pdfunite *.pdf all_plots.pdf
        '''

#-----------------------------------------------------------------------------------------
#   plot sashimi for selected SRSF genes
#-----------------------------------------------------------------------------------------


rule getIntronsForSashimi_SRSF:
    input: 
        rds = '../data/ds_v_dge/{ds_tissue_2}_v_{ds_tissue_1}_data.rds' # separate smk rules to make rds files 
    output: 'plots/sashimi/SRSF/{ds_tissue_2}_v_{ds_tissue_1}/introns.txt'
    params:
        R_script = 'workflow/scripts/getIntrons_withGene.R',
        contrast = '{ds_tissue_2}_v_{ds_tissue_1}',
        minDeltaPsi = 0,
        FDR = 1,
        minL2FC = 0,
        gene_list = 'plots/sashimi/SRSF/SRSF.genes'
    log: 'logs/getIntronsForSashimi_SRSF/{ds_tissue_2}_v_{ds_tissue_1}.log'
    shell:
        '''
        Rscript {params.R_script} {input.rds} {params.contrast} {output} {params.minDeltaPsi} {params.FDR} {params.minL2FC} {params.gene_list} &> {log}
        '''

def getPlotIntrons_SRSF(wildcards):
    with open(f'plots/sashimi/SRSF/{wildcards.ds_tissue_2}_v_{wildcards.ds_tissue_1}/introns.txt') as f:
        introns = [x.strip() for x in f.readlines()]
    return introns


use rule PrepSashimiDsGtex as PrepSashimiDsGtex_SRSF with:
    input: 
        ds_sample_group = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds_sample_group.txt',
        ds_effect_sizes = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds_effect_sizes.txt',
        ds_intron_count = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds_perind_numers.counts.noise_by_intron.gz',
        introns = 'plots/sashimi/SRSF/{ds_tissue_2}_v_{ds_tissue_1}/introns.txt',
    output:
        flag = touch('plots/sashimi/SRSF/{ds_tissue_2}_v_{ds_tissue_1}/prep.done')
    params:
        contrast = '{ds_tissue_2}_v_{ds_tissue_1}',
        outDir = 'plots/sashimi/SRSF/{ds_tissue_2}_v_{ds_tissue_1}',
        bwPrefix = 'resources/GTEx/BigWig',
        plotIntrons = getPlotIntrons_SRSF,
        iniTemplate = 'config/template-sashimi-diffsplice.ini',
        shellTemplate = 'config/template-plot-sashimi-cmd.sh',
        pyscript = 'workflow/scripts/prepSashimi.py',
    log: 'logs/PrepSashimiDsGtex/SRSF/{ds_tissue_2}_v_{ds_tissue_1}.log'


use rule PlotSashimiDsGtex as PlotSashimiDsGtex_SRSF with:
    input: 'plots/sashimi/SRSF/{ds_tissue_2}_v_{ds_tissue_1}/prep.done'
    output: touch('plots/sashimi/SRSF/{ds_tissue_2}_v_{ds_tissue_1}/plot.done')
    params:
        folder = 'plots/sashimi/SRSF/{ds_tissue_2}_v_{ds_tissue_1}',
    conda: 'pygenometracks'
    group: 'PrepSashimiDsGtex'
    threads: 4
    log: 'logs/PlotSashimiDsGtex/SRSF/{ds_tissue_2}_v_{ds_tissue_1}.log'


