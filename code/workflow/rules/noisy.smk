### Detect Noisy Introns ###



rule RunNoisy:
    input: 
        junc_list_file = 'my-temp-data/test.junclist.txt', #config['junction_list'],
        intron_class   = config['intron_class']
    output:
        perind_noise_counts = 'results/noise/out_perind.counts.noise.gz'
    params:
        run_dir    = 'results/noise',
        out_prefix = 'out', # note do not include parent dir
        py_script  = 'workflow/scripts/leafcutter_cluster_regtools_noisy_CD.py'
    threads: 1
    resources: cpu = 1, time=2000, mem=25000
    shell:
        '''
        python {params.py_script} \
            -j {input.junc_list_file} \
            -N {input.intron_class} \
            -r {params.run_dir} -o {params.out_prefix}
        '''