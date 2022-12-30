### Detect Noisy Introns ###

rule MakeIntronClustersGeuvadis:
    '''
    Use all geuvadis junctions to first construct intron cluster file
    Here, keep constitutive introns and withou any read filters
    The only filter is that introns must be within 100kb
    '''
    message: '### Make intron clusters using any GTEx samples'
    output: 
        pooled   = 'resources/{datasource}/juncs/intron_clusters/leafcutter_pooled_mkclu',
        clusters = 'resources/{datasource}/juncs/intron_clusters/leafcutter_clusters'
    params:
        run_dir    = 'resources/Geuvadis/juncs/intron_clusters',
        out_prefix = 'leafcutter',
        junc_files = '/project2/yangili1/ankeetashah/noisysQTL/geuvadis/geuvadis_juncs/ERR*Aligned.sorted.junc.hg38',
        py_script  = 'workflow/scripts/leafcutter_make_clusters.py'
    shell:
        '''
        python {params.py_script} \
            -r {params.run_dir} \
            -o {params.out_prefix} \
            -j <(realpath {params.junc_files}) \
            --includeconst
        '''


def getGeuvadisJuncs(wildcards):
    '''Get sample names ERR### based on population and available files on disk
    '''
    from glob import glob1
    ERR_Files = glob1('/project2/yangili1/ankeetashah/noisysQTL/geuvadis/geuvadis_juncs', 'ERR*Aligned.sorted.junc.hg38')
    RIDS = [e.split('.')[0] for e in ERR_Files]
    
    if wildcards.population == 'EUR': POPIDS = ['CEU', 'FIN', 'GBR', 'TSI']
    elif wildcards.population in ['CEU', 'FIN', 'GBR', 'TSI', 'YRI']: POPIDS = [wildcards.population]
    else: print(f"Error, population must be one of: 'CEU', 'FIN', 'GBR', 'TSI', 'YRI' or 'EUR'...\n")

    df = Geuvadis_Metadata.query('run_id in @RIDS and Pop_id in @POPIDS').reset_index(drop=True)
    
    sid, rid = set(list(df['Sample'])), set(list(df['run_id']))
    print(f"Getting {len(sid)} samples from {len(rid)} individuals...\n")
    print(df)
    
    return ['/project2/yangili1/ankeetashah/noisysQTL/geuvadis/geuvadis_juncs/' \
            + r + '.Aligned.sorted.junc.hg38' for r in rid]

rule AggregateJuncBySubjGeuvadis:
    '''
    Aggregate ERR based junc file into subject (person) based
    Had to name `results/pheno/aggJuncs` to differentiate from rule AggregateJuncBySubjGtex
    '''
    input:
        junc_files = getGeuvadisJuncs
    output:
        done = temp(touch('results/pheno/aggJuncs/{datasource}/{population}/done')),
        d = temp(directory('results/pheno/aggJuncs/{datasource}/{population}'))
    params:
        py_script = 'workflow/scripts/leafcutter_merge_multi_samples.py',
        clusters = 'resources/{datasource}/juncs/intron_clusters/leafcutter_clusters',
        out_prefix = 'leafcutter',
        lookup = config['Dataset']['Geuvadis']['Metadata']
    shell:
        '''
        python {params.py_script} \
            -j <(realpath {input.junc_files}) \
            -c {params.clusters} \
            -L <(awk 'BEGIN {{OFS="\t"; print "SAMPID", "SUBJID"}} NR > 1 {{print $2,$1}}' {params.lookup}) \
            -r {output.d}
        '''


rule AnnotateNoisySplicingGeuvadis:
    message:'### Annotate noisy splicing intron clusters in GTEx'
    input: 
        junc_files = 'results/pheno/aggJuncs/{datasource}/{population}',
        intron_class    = config['intron_class'],
        intron_clusters = 'resources/{datasource}/juncs/intron_clusters/leafcutter_clusters'
    output:
        perind_noise_counts = 'results/pheno/noisy-geuv/{datasource}/{population}/leafcutter_perind.counts.noise.gz'
    params:
        run_dir    = 'results/pheno/noisy-geuv/{datasource}/{population}',
        out_prefix = 'leafcutter', # note do not include parent dir
        py_script  = 'workflow/scripts/leafcutter_cluster_regtools_noisy_CD.py'
    threads: 1
    resources: cpu=1, time=2100, mem_mb=25000
    shell:
        '''
        python {params.py_script} \
            -j <(realpath {input.junc_files}/*.agg.junc.gz) \
            -r {params.run_dir} \
            -o {params.out_prefix} \
            -N {input.intron_class} \
            -c {input.intron_clusters} \
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
#     output: 'results/pheno/{datasource}/juncs/junc-file-list.txt'
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





