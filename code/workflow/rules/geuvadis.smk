### Detect Noisy Introns ###

rule MakeIntronClustersGeuvadis:
    '''
    Use all geuvadis junctions to first construct intron cluster file
    Here, keep constitutive introns and withou any read filters
    The only filter is that introns must be within 100kb
    '''
    message: '### Make intron clusters using any GTEx samples'
    output: 
        pooled   = 'resources/Geuvadis/{population}/juncs/intron_clusters/leafcutter_pooled',
        clusters = 'resources/Geuvadis/{population}/juncs/intron_clusters/leafcutter_refined_noisy'
    params:
        run_dir    = 'resources/Geuvadis/{population}/juncs/intron_clusters',
        out_prefix = 'leafcutter',
        junc_files = 'resources/Geuvadis/{population}/juncs/*.Filtered.autosomes.junc', # copied from Ben's folder
        py_script  = 'workflow/submodules/leafcutter2/scripts/leafcutter_make_clusters.py'
    shell:
        '''
        python {params.py_script} \
            -r {params.run_dir} \
            -o {params.out_prefix} \
            -j <(realpath {params.junc_files}) 
            # --includeconst
        
        ls {output.pooled}
        '''
        # commented out --includeconst on 23/6/5


rule AnnotateNoisySplicingGeuvadis:
    message:'### Annotate noisy splicing intron clusters with pre-clustering'
    input: 
        junc_files = 'resources/Geuvadis/{population}/juncs',
        intron_clusters = 'resources/Geuvadis/{population}/juncs/intron_clusters/leafcutter_refined_noisy'
    output:
        perind_noise_counts = 'results/pheno/noisy/Geuvadis/{population}/leafcutter_perind.counts.noise.gz'
    params:
        run_dir    = 'results/pheno/noisy/Geuvadis/{population}',
        out_prefix = 'leafcutter', # file name prefix, not directory
        intron_class = ','.join(config['intron_class']),
        py_script  = 'workflow/submodules/leafcutter2/scripts/leafcutter2_regtools.py'
    threads: 1
    resources: cpu=1, time=2100, mem_mb=25000
    shell:
        '''
        python {params.py_script} \
            -j <(realpath {input.junc_files}/*.junc) \
            -r {params.run_dir} \
            -o {params.out_prefix} \
            -N {params.intron_class} \
            -c {input.intron_clusters} \
            -k 
        '''


rule AnnotateNoisySplicingGeuvadis_Wo_PreCluster:
    message: '### Annotate noisy splicing intron clusters without pre-clustering'
    input:
        junc_files = 'resources/Geuvadis/{population}/juncs'
    output:
        perind_noise_counts = 'results/pheno/noisy/Geuvadis/{population}/leafcutter_wo_precluster_perind.counts.noise.gz'
    params:
        run_dir    = 'results/pheno/noisy/Geuvadis/{population}',
        out_prefix = 'leafcutter_wo_precluster', # file name prefix, not directory,
        intron_class = ','.join(config['intron_class']),
        py_script  = 'workflow/submodules/leafcutter2/scripts/leafcutter_cluster_regtools_noisy_CD_v2.py'
    threads: 1
    resources: cpu=1, time=2100, mem_mb=25000
    shell:
        '''
        python {params.py_script} \
            -j <(realpath {input.junc_files}/*.junc) \
            -r {params.run_dir} \
            -o {params.out_prefix} \
            -N {params.intron_class} \
            -k 
        '''





# raw junctions from ankeeta's folder
# /project2/yangili1/ankeetashah/noisysQTL/geuvadis/geuvadis_juncs/ERR188021.Aligned.sorted.junc.hg38

# rule MakeJunctionList:
#     '''
#     Procedures:
#         1.  get all junction files from directory
#         2.  subset samples that have genotype and only have 1 ERR ID per sample
#         3.  construct junction file path of such samples
#         4.  save file paths in a file
#     '''
#     message: '### Make a list of junction input file'
#     output: 'results/pheno/Geuvadis/juncs/junc-file-list.txt'
#     params:
#         junc_dir = '/project2/yangili1/ankeetashah/noisysQTL/geuvadis/geuvadis_juncs/'
#     run:
#         # junction file names
#         junc_files = glob.glob1(params.junc_dir, 'ERR*.Aligned.sorted.junc.hg38')
#         junc_ERRs= [j.split(".")[0] for j in junc_files] # ERR IDs in junc files
        
#         # junction ERRs corresponding Sample IDs
#         junc_samples = list(set(Geuvadis_Metadata.query('run_id in @junc_ERRs').Sample))

#         # intersect junc file representing sample IDs with 1to1 genotyped sample IDs
#         junc_keep_samples = [x for x in junc_samples if x in Geuvadis_Linked_SampleIDs]

#         # ERR IDs that are in genotyped 1to1 samples and junction file samples
#         junc_keep_ERRs = list(Geuvadis_Metadata.query('Sample in @junc_keep_samples').run_id)
#         junc_keep_ERRs = [c for c in junc_ERRs if c in junc_keep_ERRs ]

#         # full file path of selected junction files
#         junc_file_paths = [params.junc_dir + c + '.Aligned.sorted.junc.hg38' for c in junc_keep_ERRs]

#         # write out paths in a file
#         with open(output[0], 'w') as f:
#             f.write('\n'.join(junc_file_paths))





