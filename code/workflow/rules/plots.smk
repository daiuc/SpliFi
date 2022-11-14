

rule Sushi:
    message: 'plot sashimi plots'
    output: 'mysashimi.pdf'
    singularity: '/scratch/midway2/chaodai/singularity/ggsashimi.sif'
    params:
        bams = 'workflow/submodules/ggsashimi/examples/input_bams.tsv',
        region = 'chr10:27040584-27048100',
        gtf = '/project2/yangili1/cdai/genome_index/hs38/gencode.v38.primary_assembly.annotation.gtf.gz',
        outprefix = 'mysashimi',
        pyscript = 'workflow/submodules/ggsashimi/ggsashimi.py'
    shell:
        '''
        {params.pyscript} -b {params.bams} -c {params.region} -g {params.gtf} -o {params.outprefix}
        ls -l {output}
        '''

