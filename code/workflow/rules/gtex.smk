### NOTE:
#   -   the `tissue` wildcards have `-` replacing white spaces. While the `SMTS`
#       and the `SMTSD` fields has white space, that's why in rules using `tissue`
#       wildcard, there must be a step to replace `-` back to space.


rule MungeGTExJuncs:
    message: '# Munge GTEx junc file'
    input: 'resources/GTEx/juncs/GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct.gz'
    output: 
        'resources/GTEx/juncs/munged_juncs.txt.gz'
    params: 
        anno = config['annotation']['gencode_v26_genes'],
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

# since all GTEx sample files have the same introns, only need to run this once
rule MakeIntronClustersGtex:
    message: '### Make intron clusters using any GTEx samples'
    output: 
        pooled   = 'resources/GTEx/juncs/intron_clusters/GTEx_pooled_mkclu',
        clusters = 'resources/GTEx/juncs/intron_clusters/GTEx_clusters'
    params:
        run_dir    = 'resources/GTEx/juncs/intron_clusters',
        out_prefix = 'GTEx',
        junc_files = 'resources/GTEx/juncs/split_juncs',
        py_script  = 'workflow/scripts/leafcutter_make_clusters.py'
    shell:
        '''
        python {params.py_script} -r {params.run_dir} -o {params.out_prefix} \
            -j <(realpath {params.junc_files}/*.tsv.gz | head)
        '''


# aggregate junction file by subject id
rule AggregateJuncBySubjGtex:
    output: 
        out_flag = temp(touch('results/pheno/groupJuncs/{datasource}/{tissue}/groupjunc.done')),
        out_dir = temp(directory('results/pheno/groupJuncs/{datasource}/{tissue}/'))
    params:
        junc_prefix = 'resources/GTEx/juncs/split_juncs/',
        junc_suffix = '.tsv.gz',
        out_prefix  = 'results/pheno/groupJuncs/{datasource}/{tissue}/',
        out_suffix  = '.tsv'
    run:
        from functools import reduce

        tissue = wildcards.tissue.replace("-", " ")
        if tissue in list(Gtex_Metadata.SMTS):
            df = Gtex_Metadata.query('SMTS in @tissue').reset_index(drop=True).drop_duplicates()
        elif tissue in list(Gtex_Metadata.SMTSD):
            df = Gtex_Metadata.query('SMTSD in @tissue').reset_index(drop=True).drop_duplicates()
        else:
            print("Error. Check the tissue type you entered!")
            exit(1)
        samps = {} # {subjid: [sampid, sampid, ...]}
        for k,v in df.groupby('SUBJID'):
            samps[k] = list(v.SAMPID)

        print(f'Make junction files for Tissue type: {tissue}')
        print(f'Total {len(samps.keys())} individuals, {len([x for l in list(samps.values()) for x in l])} samples ...\n')

        for k, v in samps.items():
            junc_files = [params.junc_prefix + x + params.junc_suffix for x in v] # junc files
            out_file = output['out_dir'] + '/' + k + params.out_suffix
            cts = [] # list of list for counts
            for j in junc_files:
                with gzip.open(j) as f:
                    lines = [ln.decode().strip().split() for ln in f.readlines()]
                    ct = [int(x[4]) for x in lines]
                    cts.append(ct)
            cts = reduce(lambda x,y: [a+b for a,b in zip(x,y)], cts) # sum reads by SUBJID if having multiple SAMPID
            buf = ['\t'.join([lines[i][0], lines[i][1], lines[i][2], 
                              lines[i][3], str(cts[i]), lines[i][5]]
                            ) + '\n' for i in range(len(lines))]
            with open(out_file, 'w') as outf:
                outf.writelines(buf)
                # print(f'Wring {k} ... Done.')
                
            

rule AnnotateNoisySplicingGtex:
    message:'### Annotate noisy splicing intron clusters in GTEx'
    input: 
        junc_file_dir  = 'results/pheno/groupJuncs/{datasource}/{tissue}',
        intron_class    = config['intron_class'],
        intron_clusters = 'resources/GTEx/juncs/intron_clusters/GTEx_clusters'
    output:
        perind_noise_counts = 'results/pheno/noisy/{datasource}/{tissue}/leafcutter_perind.counts.noise.gz'
    params:
        run_dir    = 'results/pheno/noisy/{datasource}/{tissue}',
        out_prefix = 'leafcutter', # note do not include parent dir
        py_script  = 'workflow/scripts/leafcutter_cluster_regtools_noisy_CD.py'
    threads: 1
    resources: cpu=1, time=2100, mem_mb=25000
    shell:
        '''
        python {params.py_script} \
            -j <(realpath {input.junc_file_dir}/*.tsv) \
            -r {params.run_dir} \
            -o {params.out_prefix} \
            -N {input.intron_class} \
            -c {input.intron_clusters} \
            -k 
        '''
