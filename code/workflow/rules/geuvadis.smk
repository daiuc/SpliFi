'''
NOTE:
    - As of 2024/1/15, only used the EUR populations, which include CEU, FIN, GBR, TSI etc.
    - junction files are copied from Ben's project folder /project2/yangili1/bjf79/ChromatinSplicingQTLs/code/Alignments/STAR_Align/Expression.Splicing
      I used `code/resources/get-geuvadis-juncs.R.ipynb` to copy the junction files to `resources/Geuvadis/{population}/juncs`
    
'''


# optional rule used to speed up the workflow, intron clusters may run once
# and reuse as needed.
rule ClusterIntronGeuvadis:
    message: '### Make intron clusters'
    output: 
        pooled   = 'results/pheno/noisy/Geuvadis/{population}/leafcutter_pooled',
        clusters = 'results/pheno/noisy/Geuvadis/{population}/leafcutter_refined_noisy',
        lowusage = 'results/pheno/noisy/Geuvadis/{population}/leafcutter_lowusage_introns', # intermediate
        refined  = 'results/pheno/noisy/Geuvadis/{population}/leafcutter_refined' # intermediate
    params:
        run_dir    = 'results/pheno/noisy/Geuvadis/{population}',
        out_prefix = 'leafcutter',
        junc_files = 'resources/Geuvadis/{population}/juncs/*.Filtered.autosomes.junc', # copied from Ben's folder
        py_script  = 'workflow/submodules/leafcutter2/scripts/leafcutter_make_clusters.py'
    log: 'logs/ClusterIntronGeuvadis/{population}.log'
    shell:
        '''
        python {params.py_script} \
            -r {params.run_dir} \
            -o {params.out_prefix} \
            -j <(realpath {params.junc_files}) &> {log}
        ls -lah {output.pooled} {output.clusters} {output.lowusage} {output.refined} &>> {log}

        '''


rule Leafcutter2Geuvadis:
    message:'### Run leafcutter2 on Geuvadis samples'
    input: 
        junc_files = 'resources/Geuvadis/{population}/juncs',
        intron_clusters = 'results/pheno/noisy/Geuvadis/{population}/leafcutter_refined_noisy',
    output:
        perind_noise_by_intron = 'results/pheno/noisy/Geuvadis/{population}/leafcutter_perind.counts.noise_by_intron.gz'
    params:
        py_script  = 'workflow/submodules/leafcutter2/scripts/leafcutter2_regtools.py',
        run_dir    = 'results/pheno/noisy/Geuvadis/{population}',
        out_prefix = 'leafcutter', # note do not include parent dir
        gtf = config['annotation']['gtf']['v43'],
        genome = config['genome38'],
        pre_clustered = '-c results/pheno/noisy/Geuvadis/{population}/leafcutter_refined_noisy', 
        other_params = '-k' # not keeping constitutive introns
    threads: 1
    resources: cpu=1, time=2100, mem_mb=25000
    log: 'logs/Leafcutter2Geuvadis/{population}.log'
    shell:
        '''
        python {params.py_script} \
            -j <(realpath {input.junc_files}/*.junc) \
            -r {params.run_dir} \
            -o {params.out_prefix} \
            -A {params.gtf} \
            -G {params.genome} \
            {params.pre_clustered} {params.other_params} &> {log}
        
        ls {output.perind_noise_by_intron} &>> {log}

        '''



use rule Leafcutter2Geuvadis as Leafcutter2Geuvadis_wConst with:
    message: '### run leafcutter2 on Geuvadis samples with constitutive introns'
    input: 
        junc_files = 'resources/Geuvadis/{population}/juncs'
    output:
        perind_noise_by_intron = 'results/pheno/noisy/Geuvadis/{population}/wConst_perind.constcounts.noise_by_intron.gz'
    params:
        py_script  = 'workflow/submodules/leafcutter2/scripts/leafcutter2_regtools.py',
        run_dir    = 'results/pheno/noisy/Geuvadis/{population}',
        pre_clustered = '',
        out_prefix = 'wConst',
        gtf = config['annotation']['gtf']['v43'],
        genome = config['genome38'],
        other_params = '-k --includeconst' # include constitutive introns
    log: 'logs/Leafutter2Geuvadis_wConst/{population}.log'



rule Leaf2AnnotIntrons:
    '''
    Script is modified from leafcutter2_regtools.py to just use the annotation step
    to annotate all intron junctions in dataset.
    Used for paper Figure1C
    '''
    message: '### Use leafcutter2 to annotate ALL introns (no filtering)'
    input:
        perind_counts = 'results/pheno/noisy/Geuvadis/EUR/wConst_perind.constcounts.gz',
        junc_files = 'results/pheno/noisy/Geuvadis/EUR/wConst_junction_classifications.txt',
    output:
        noisydiag = f"results/pheno/noisy/Geuvadis/EUR/wConst_perind.constcounts.annotated.gz",
        numersdiag = f"results/pheno/noisy/Geuvadis/EUR/wConst_perind.numers.annotated.gz"
    params:
        out_prefix = 'results/pheno/noisy/Geuvadis/EUR/wConst_perind',
        py_script = 'workflow/scripts/lf2_label_introns.py'
    shell:
        '''
        python {params.py_script} \
            -i {input.perind_counts} \
            -o {params.out_prefix} \
            -a {input.junc_files} \

        ls {output}

        '''


# rule AnnotateNoisySplicingGeuvadis_const:
#     input:
#         junc_files = 'resources/Geuvadis/{population}/juncs'
#     output:
#         perind_noise_counts = 'results/pheno/noisy/Geuvadis/{population}/allIntron/allIntron_perind.constcounts.noise.gz'
#     params:
#         run_dir    = 'results/pheno/noisy/Geuvadis/{population}/allIntron',
#         out_prefix = 'allIntron', # file name prefix, not directory,
#         intron_class = ','.join(config['intron_class']),
#         py_script  = 'workflow/submodules/leafcutter2/scripts/leafcutter2_regtools.py'
#     threads: 1
#     resources: cpu=1, time=2100, mem_mb=25000
#     shell:
#         '''
#         python {params.py_script} \
#             -j <(realpath {input.junc_files}/*.junc) \
#             -r {params.run_dir} \
#             -o {params.out_prefix} \
#             -N {params.intron_class} \
#             -k --includeconst
#         '''






