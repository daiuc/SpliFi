### Detect Noisy Introns ###


# raw junctions from ankeeta's folder
# /project2/yangili1/ankeetashah/noisysQTL/geuvadis/geuvadis_juncs/ERR188021.Aligned.sorted.junc.hg38

rule MakeJunctionList:
    '''
    Procedures:
        1.  get all junction files from directory
        2.  subset samples that have genotype and only have 1 ERR ID per sample
        3.  construct junction file path of such samples
        4.  save file paths in a file
    '''
    message: '### Make a list of junction input file'
    output: 'results/pheno/{datasource}/juncs/junc-file-list.txt'
    params:
        junc_dir = '/project2/yangili1/ankeetashah/noisysQTL/geuvadis/geuvadis_juncs/'
    run:
        # junction file names
        junc_files = glob.glob1(params.junc_dir, 'ERR*.Aligned.sorted.junc.hg38')
        junc_ERRs= [j.split(".")[0] for j in junc_files] # ERR IDs in junc files
        
        # junction ERRs corresponding Sample IDs
        junc_samples = list(set(Geuvadis_Metadata.query('run_id in @junc_ERRs').Sample))

        # intersect junc file representing sample IDs with 1to1 genotyped sample IDs
        junc_keep_samples = [x for x in junc_samples if x in Geuvadis_Linked_SampleIDs]

        # ERR IDs that are in genotyped 1to1 samples and junction file samples
        junc_keep_ERRs = list(Geuvadis_Metadata.query('Sample in @junc_keep_samples').run_id)
        junc_keep_ERRs = [c for c in junc_ERRs if c in junc_keep_ERRs ]

        # full file path of selected junction files
        junc_file_paths = [params.junc_dir + c + '.Aligned.sorted.junc.hg38' for c in junc_keep_ERRs]

        # write out paths in a file
        with open(output[0], 'w') as f:
            f.write('\n'.join(junc_file_paths))


rule RunNoisy:
    message:'### Detect noisy splicing.'
    input: 
        junc_list_file =  rules.MakeJunctionList.output,
        intron_class   = config['intron_class']
    output:
        perind_noise_counts = 'results/pheno/{datasource}/noise_counts/{datasource}_perind.counts.noise.gz'
    params:
        run_dir    = 'results/pheno/{datasource}/noise_counts',
        out_prefix = '{datasource}', # note do not include parent dir
        py_script  = 'workflow/scripts/leafcutter_cluster_regtools_noisy_CD.py'
    threads: 1
    resources: cpu=1, time=2100, mem_mb=25000
    shell:
        '''
        python {params.py_script} \
            -j {input.junc_list_file} \
            -N {input.intron_class} \
            -r {params.run_dir} -o {params.out_prefix}
        '''