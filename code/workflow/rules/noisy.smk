### Detect Noisy Introns ###


# process GTEx





def getRunNoisyInputs(wildcards):
    return config['junction_list'].get(wildcards.datasource)

rule RunNoisy:
    input: 
        junc_list_file =  getRunNoisyInputs,
        intron_class   = config['intron_class']
    output:
        perind_noise_counts = 'results/{datasource}/noise/{datasource}_perind.counts.noise.gz'
    wildcard_constraints:
        datasource = '(test)|(Geuvadis)|(GTEx)|(TCGA)'
    params:
        run_dir    = 'results/{datasource}/noise',
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