

rule FeatureCountGeneExpression:
    message: 'FeatureCount Gene expression'
    input: 
        anno = config['annotation']['gencode_v38_genes'],
        bam = '/project2/yangili1/ankeetashah/hg38_LCL/{geuvadis_sid}.bam'
    output: 
        summary = 'results/pheno/{datasource}/{population}/ge/{geuvadis_sid}.counts.summary'
    params: 
        TMP = '/home/chaodai/scratch/TMP',
        prefix = 'results/pheno/{datasource}/{population}/ge/{geuvadis_sid}'
    threads: 4
    resources: cpu = 4, mem_mb = 15000, time = 2100
    shell: 
        '''
        featureCounts \
            -F SAF --minOverlap 10 \
            -p --countReadPairs -P -B -d 50 -D 20000 \
            -T {threads} --tmpDir {params.TMP} \
            -Q 20 -M \
            -a <(awk 'BEGIN {{OFS="\t"}}; NR>1 {{print $4,$1,$2,$3,$6}}' {input.anno}) \
            -o {params.prefix}.counts {input.bam}
        '''