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
        with open('resources/GTEx/juncs/SampleIDs.txt') as f:
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
