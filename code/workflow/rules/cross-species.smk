'''
Data downloaded from here: 
- Raw and processed RNA-seq data have been deposited in ArrayExpress with the accession codes: E-MTAB-6769 (chicken), E-MTAB-6782 (rabbit), E-MTAB-6798 (mouse), E-MTAB-6811 (rat), E-MTAB-6813 (rhesus macaque), E-MTAB-6814 (human) and E-MTAB-6833 (opossum) (https://www.ebi.ac.uk/arrayexpress/). The temporal profiles of individual genes across organs and species can be visualized and downloaded using the web-based application: http://evodevoapp.kaessmannlab.org.

Dowloaded 10 bams from each species and tissue using globus. Download shell script within each folder: eg. `resources/moreira/bam/Human/globus.sh`


Species	Genome assembly	Annotation used in GSNAP	Annotation used by HTSeq to generate counts per gene
Gallus gallus	Galgal4	71	77
Mus musculus	GRCm38 (mm10)	71	77
Homo sapiens	GRCh37 (hg19)	71	75
Monodelphis domestica	BROADO5	71	77
Rattus norvegicus	Rnor_5.0	73	77
Oryctolagus cuniculus	OryCun2.0	75	77
Rhesus macaque	MMUL_1	77	77
'''


import pandas as pd
import glob

configfile: "config/config.yaml"

AUTOSOMES = ['chr'+str(i) for i in range(1,23)]
CHROMS = AUTOSOMES + ['chrX', 'chrY']


# remove white space and fix parentethies in tissue name
fixChr = str.maketrans({" ": None, "(": "_", ")": "_"})


wildcard_constraints:
    specie = "Human|Mouse|Rat|Macaque",
    tissue = "Brain|Liver|Heart|Kidney|Testis",


allf = glob.glob('resources/moreira/bam/*/*bam')
allf = [x.split('/')[-1].split('.') for x in allf]
allf = [(x[1], x[2], x[0]) for x in allf]
sampleInfo = pd.DataFrame(allf, columns=['specie', 'tissue', 'sample']).sort_values(['specie', 'tissue', 'sample']).drop_duplicates()


contrasts = ["Liver_v_Kindey", "Brain_v_Liver", "Liver_v_Testis", "Heart_v_Liver", "Brain_v_Kidney", "Kidney_v_Testis", "Heart_v_Kidney", "Brain_v_Testis", "Brain_v_Heart", "Heart_v_Testis"]


localrules: all

rule all:
    input:
        expand("resources/moreira/juncs/{specie}/{tissue}-{sample}.junc", zip, 
               specie=sampleInfo['specie'].to_list(),
               tissue=sampleInfo['tissue'].to_list(),
               sample=sampleInfo['sample'].to_list()
               ),
        # expand("results/ds/multiSpecies/{specie}/{contrasts}/ds_sample_group.txt",
        #        specie=sampleInfo['specie'].unique(),
        #        contrasts=contrasts
        #        ),

def getBam(wildcards):
    specie, tissue, sample = wildcards.specie, wildcards.tissue, wildcards.sample
    bams = glob.glob(f'resources/moreira/bam/{specie}/*{specie}.{tissue}*bam')
    bams = [b for b in bams if specie in b and tissue in b]
    s = [x.split('/')[-1].split('.')[0] for x in bams]
    lookup = {k: v for k, v in zip(s, bams)}

    return lookup[sample]

rule SortIndexBam: 
    input:
        bam = getBam
    output: 
        bam = temp("resources/moreira/sortedbam/{specie}/{tissue}-{sample}.bam"),
        bai = temp("resources/moreira/sortedbam/{specie}/{tissue}-{sample}.bam.bai")
    params:
        chroms = lambda x: ','.join(CHROMS)
    threads: 4
    group: "bam-junc"
    shell:
        '''
        samtools sort -@ {threads} -o {output.bam} {input.bam}
        samtools index -@ {threads} {output.bam}
        '''

rule ExtractIntrons:
    input: 
        bam = "resources/moreira/sortedbam/{specie}/{tissue}-{sample}.bam",
        bai = temp("resources/moreira/sortedbam/{specie}/{tissue}-{sample}.bam.bai"),
    output:
        'resources/moreira/juncs/{specie}/{tissue}-{sample}.junc'
    threads: 1
    resources:
        cpu = 1, mem_mb = 25000, time = 2000
    group: "bam-junc"
    shell:
        '''
        regtools junctions extract -s XS -o {output} {input.bam}
        ls {output}
        '''

# ----------------------------------------------------------------------------------------
#           Cluster introns using all GTEx tissues at once
#         Then use the clusters to run leafcutter2 on specific tissues
#         This enables differential splicing analysis across tissues
# ----------------------------------------------------------------------------------------


# optional rule used to speed up the workflow, intron clusters may run once
# and reuse as needed.
rule ClusterIntrons:
    output: 
        pooled   = 'results/ds/multiSpecies/{specie}/all_tissues_pooled',
        clusters = 'results/ds/multiSpecies/{specie}/all_tissues_refined_noisy',
        lowusage = 'results/ds/multiSpecies/{specie}/all_tissues_lowusage_introns', # intermediate
        refined  = 'results/ds/multiSpecies/{specie}/all_tissues_refined' # intermediate
    params:
        run_dir    = 'results/ds/multiSpecies/{specie}',
        out_prefix = 'all_tissues',
        junc_files = 'resources/moreira/juncs/{specie}',
        py_script  = 'workflow/submodules/leafcutter2/scripts/leafcutter_make_clusters.py'
    log:'logs/crossSpecies/ClusterIntrons{specie}.log'
    resources: cpu=1, time=2100, mem_mb=35000
    shell:
        '''
        python {params.py_script} \
            -r {params.run_dir} \
            -o {params.out_prefix} \
            -j <(ls {params.junc_files}/*.junc) &> {log}
        
        ls -lah {output.pooled} {output.clusters} {output.lowusage} {output.refined} &>> {log}
        
        '''


rule LeafcutterForDSGtex:
    message: '### Run leafcutter2 on GTEx samples for differential splicing analysis'
    input:
        pre_clusters = 'results/ds/multiSpecies/{specie}/all_tissues_refined_noisy',
    output:
        ds_numers = 'results/ds/multiSpecies/{specie}/{ds_tissue_2}_v_{ds_tissue_1}/ds_perind_numers.counts.noise_by_intron.gz',
        # ds_counts_lf1 is necessary for leafcutter_ds.R script
        ds_numers_lf1 = 'results/ds/multiSpecies/{specie}/{ds_tissue_2}_v_{ds_tissue_1}/ds_perind_numers.counts.noise_by_intron.lf1.gz',
        ds_sample_group = 'results/ds/multiSpecies/{specie}/{ds_tissue_2}_v_{ds_tissue_1}/ds_sample_group.txt'
    params:
        run_dir = 'results/ds/multiSpecies/{specie}/{ds_tissue_2}_v_{ds_tissue_1}',
        out_prefix = 'ds',
        tissue1_juncs = 'resources/moreira/juncs/{specie}/{ds_tissue_1}',
        tissue2_juncs = 'resources/moreira/juncs/{specie}/{ds_tissue_2}',
        NSamples = 10, # select this number of samples per tissue type
        pre_clustered = '-c results/ds/multiSpecies/{specie}/all_tissues_refined_noisy',
        gtf = lambda w: config['crossspecies']['gtf'][w.specie], 
        genome = lambda w: config['crossspecies']['genome'][w.specie],
        other_params = '-k ', # not keeping constitutive introns
        py_script  = 'workflow/submodules/leafcutter2/scripts/leafcutter2_regtools.py',
        py_script2 = 'workflow/scripts/makeSampleGroupFileForDifferentialSplicing.py'
    log: 'logs/cross-species/LeafcutterForDSGtex/{specie}/{ds_tissue_2}_v_{ds_tissue_1}.log'
    resources: cpu = 1, time = 2100, mem_mb = 25000
    shell:
        '''
        python {params.py_script} \
            -j <(cat <(ls {params.tissue1_juncs}*.junc ) <(ls {params.tissue2_juncs}*.junc)) \
            -r {params.run_dir} \
            -o {params.out_prefix} \
            -A {params.gtf} \
            -G {params.genome} \
            {params.pre_clustered} {params.other_params} &> {log}

        # make sample group file for differential splicing analysis
        python {params.py_script2} -i {output.ds_numers} -o {output.ds_numers_lf1} -s {output.ds_sample_group} &>> {log}

        '''


# rule RunLeafcutterDiffSplicingGtex:
#     input:
#         ds_numers_lf1 = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds_perind_numers.counts.noise_by_intron.lf1.gz',
#         ds_sample_group = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds_sample_group.txt'
#     output:
#         flag = touch('results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/done')
#         # produces two files:
#         # 1. {outprefix}_effect_sizes.txt
#         # 2. {outprefix}_manual_ds_cluster_significance.txt
#     params:
#         Rscript = 'workflow/submodules/leafcutter/scripts/leafcutter_ds.R', 
#         outprefix = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds', # note you need to include path!
#         MIN_SAMPLES_PER_INTRON = 5,
#         MIN_SAMPLES_PER_GROUP = 3,
#         MIN_COVERAGE = 5
#     resources: cpu = 4, mem_mb = 30000, time = 2100
#     threads: 4
#     conda: 'leafcutter1'
#     log: 'logs/RunLeafcutterDiffSplicingGtex/{ds_tissue_2}_v_{ds_tissue_1}.log'
#     shell:
#         '''
#         {params.Rscript} --num_threads {threads} \
#             --output_prefix {params.outprefix} \
#             --min_samples_per_intron={params.MIN_SAMPLES_PER_INTRON} \
#             --min_samples_per_group={params.MIN_SAMPLES_PER_GROUP} \
#             --min_coverage={params.MIN_COVERAGE} \
#             {input.ds_numers_lf1} {input.ds_sample_group} &> {log}
#
#
#         '''
#
# # NOTE to make leafcutter's leafcutter_ds.R script work, 
# # must modify the _perind_numers.counts.noise_by_intron.gz file
# # First row should not have the 'chrom' first column
# # first columns can only be like `chr1:827775:829002:clu_1_+` because 
# # the R function only expect to split by ":" into 4 columns.
#
#
#
# # ----------------------------------------------------------------------------------------
# # ad hoc
#
# rule adhoc_test_ds_step1: # testing differential splicing with neg control
#     input: 
#         ds_numers_lf1 = 'results/ds/GTEx/Brain-Cerebellum_v_Liver/ds_perind_numers.counts.noise_by_intron.lf1.gz',
#     output:
#         neg_control_numers = 'results/ds/GTEx/ds_test/BC_v_Liver/perind.numers.gz',
#         neg_control_groups = 'results/ds/GTEx/ds_test/BC_v_Liver/sample_group.txt'
#     run:
#         import gzip
#         outf1 = gzip.open(output.neg_control_numers, 'wt')
#         outf2 = open(output.neg_control_groups, 'w')
#         with gzip.open(input.ds_numers_lf1, 'rt') as f:
#             i = 0
#             for ln in f:
#                 if i == 0:
#                     header = ln.split()
#                     cols = ([x.replace('Brain-Cerebellum', 'BC-group1') for x in header[0:100]] +
#                             [x.replace('Brain-Cerebellum', 'BC-group2') for x in header[100:200]])
#                     groups = ['BC-group1' for x in header[0:100]] + ['BC-group2' for x in header[100:200]]
#                     for c,g in zip(cols, groups):
#                         outf2.write(f'{c} {g}\n')
#                     outf1.write(' '.join(cols) + '\n')
#                 if i > 0:
#                     outln = ln.split()[:201]
#                     outf1.write(' '.join(outln) + '\n')
#                 i += 1
#
#         outf1.close()
#         outf2.close()
#
#
# use rule RunLeafcutterDiffSplicingGtex as adhoc_test_ds_step2 with:
#     input:
#         ds_numers_lf1 = 'results/ds/GTEx/ds_test/BC_v_Liver/perind.numers.gz',
#         ds_sample_group = 'results/ds/GTEx/ds_test/BC_v_Liver/sample_group.txt'
#     output:
#         flag = touch('results/ds/GTEx/ds_test/BC_v_Liver/ds.done')
#         # produces two files:
#         # 1. {outprefix}_effect_sizes.txt
#         # 2. {outprefix}_manual_ds_cluster_significance.txt
#     params:
#         Rscript = 'workflow/submodules/leafcutter/scripts/leafcutter_ds.R', 
#         outprefix = 'results/ds/GTEx/ds_test/BC_v_Liver/ds', # note you need to include path!
#         MIN_SAMPLES_PER_INTRON = 5,
#         MIN_SAMPLES_PER_GROUP = 3,
#         MIN_COVERAGE = 5
#     log: 'results/ds/GTEx/ds_test/BC_v_Liver/log'
#
#
#
# ## -----------------------------------------------------------------------------
# ##   GTEx expression data
# ## -----------------------------------------------------------------------------
#
# rule ExtractGTExGeneExpression:
#     input: 
#       tpm = 'resources/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz',
#       cnt = 'resources/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz'
#     output: 
#       tpm = 'resources/GTEx/expression/{tissue}_gene_tpm.tsv.gz',
#       cnt = 'resources/GTEx/expression/{tissue}_gene_reads.tsv.gz',
#     params: 
#         py_script = 'workflow/scripts/extract_gtex_gene_expression.py',
#         junc_meta = config['Dataset']['GTEx']['Junc_meta']
#     log: 'logs/ExtractGTExGeneExpression/{tissue}.log'
#     shell:
#         '''
#         echo extracting tpm from {input.tpm} ...&> {log}
#         python {params.py_script} \
#             -I {input.tpm} \
#             -M {params.junc_meta} \
#             -O {output.tpm} \
#             -T {wildcards.tissue} &>> {log}
#
#         echo extracting raw counts from {input.cnt} ...
#         python {params.py_script} \
#             -I {input.cnt} \
#             -M {params.junc_meta} \
#             -O {output.cnt} \
#             -T {wildcards.tissue} &>> {log}
#
#         echo "Number of lines in {output.tpm} : $(zcat {output.tpm} | wc -l)" &>> {log}
#         echo "Number of lines in {output.cnt} : $(zcat {output.cnt} | wc -l)" &>> {log}
#         '''
#
# # prepare GTEX dge data
# # NOTE: need to rerun this step and below because right now samples
# # are not matched with ds samples! 
# # Should use the ds_sample_group.txt file to get matching samples for each tissue
# rule PrepareGTExDGE:
#     input: 
#         cnt1 = 'resources/GTEx/expression/{dge_tissue1}_gene_reads.tsv.gz',
#         cnt2 = 'resources/GTEx/expression/{dge_tissue2}_gene_reads.tsv.gz',
#         ds_samples = 'results/ds/GTEx/{dge_tissue2}_v_{dge_tissue1}/ds_sample_group.txt'
#     output: 
#         #NOTE: tissue2 is intended to be numerator, and tissue 1 denominator, in subsequent dge step
#         cnt = 'results/dge/GTEx/{dge_tissue2}_v_{dge_tissue1}_counts.tsv',
#         coldata = 'results/dge/GTEx/{dge_tissue2}_v_{dge_tissue1}_coldata.tsv',
#     params:
#         R_script = 'workflow/scripts/prepare_GTEx_dge.R',
#         outdir = 'results/dge/GTEx'
#     log: 'logs/PrepareGTExDGE/{dge_tissue2}_v_{dge_tissue1}.log'
#     shell:
#         '''
#         Rscript {params.R_script} {input.cnt1} {input.cnt2} {input.ds_samples} {params.outdir} &> {log}
#         ls {output.cnt} {output.coldata} &>> {log}
#         '''
#
# # run dge
# rule DgeGtex:
#   input: 
#     cnt = 'results/dge/GTEx/{dge_tissue2}_v_{dge_tissue1}_counts.tsv',
#     coldata = 'results/dge/GTEx/{dge_tissue2}_v_{dge_tissue1}_coldata.tsv'
#   output:
#     dge = 'results/dge/GTEx/{dge_tissue2}_v_{dge_tissue1}_dge_genes.tsv',
#   params:
#     R_script = 'workflow/scripts/dge.R',
#     outprefix = 'results/dge/GTEx/{dge_tissue2}_v_{dge_tissue1}',
#     min_reads = 10,
#     min_samples = 10,
#   log: 'logs/DgeGtex/{dge_tissue2}_v_{dge_tissue1}.log'
#   shell:
#     '''
#     Rscript {params.R_script} \
#             {input.cnt} {input.coldata} {params.outprefix} \
#             {params.min_reads} {params.min_samples}  &> {log}
#     ls {output.dge} &>> {log}
#
#     '''




