'''
Prepare data used for notebooks in analysis folder. In general, the purpose of
this script is to prepare data that requires a long time to compute. Thus I 
preprocess the data and save it to a file (often a rds file).
'''

import glob
import pandas as pd



# Get list of tissues with at least 80 samples
gtex_meta = pd.read_csv('code/resources/GTEx/juncs/sampid-smts-smtsd-subjid.tsv', sep='\t')

# remove white space and fix parentethies in tissue name
fixChr = str.maketrans({" ": None, "(": "_", ")": "_"})

chosen_tissues = (
    gtex_meta.groupby(["SMTS", "SMTSD"])["SUBJID"]
    .nunique()
    .sort_values(ascending=False)
    .reset_index()
    .query("SUBJID >= 80")['SMTSD']
    .unique().tolist()
)

chosen_tissues = [x.translate(fixChr) for x in chosen_tissues]

# tissue pairs for differential analysis
with open('code/config/ds_dge_tissue_contrast_levels.txt') as f:
    CONTRASTS = f.read().splitlines()


localrules:
    GetEqtlDataGeuvadis,


rule all:
    input: 
        # extract numerators from leafcutter2 with const run's results
        #expand('data/ExtractFractions/GTEx/{tissue}.numerators_constcounts.noise_by_intron.rds',
        #       tissue=chosen_tissues),
        # categorize sQTLs into u-sQTLs and p-sQTLs, results in tsv file
        # expand('data/CategorizeSplicingQTLs/GTEx/{tissue}.tsv',
        #        tissue=chosen_tissues),
        # get all dge and ds results and save to rds files
        expand('data/ds_v_dge/{a_v_b}_data.rds', a_v_b=CONTRASTS),


rule ExtractNumFromFracGTEx:
    message: 'Extract numerator from junction reads from leafcutter2 with const counts'
    input: 'code/results/pheno/noisy/GTEx/{tissue}/wConst_perind.constcounts.noise_by_intron.gz'
    output: 'data/ExtractFractions/GTEx/{tissue}.numerators_constcounts.noise_by_intron.rds'
    params:
        script = 'scripts/getleafcuttercountmatrix2.R',
    shell:
        '''
        Rscript {params.script} {input} {output}
        '''

use rule ExtractNumFromFracGTEx as ExtractNumFromFracGeuvadis with:
    input: 'code/results/pheno/noisy/Geuvadis/{population}/wConst_perind.constcounts.noise_by_intron.gz'
    output: 'data/ExtractFractions/Geuvadis/{population}.numerators_constcounts.noise_by_intron.rds'


rule CategorizeSplicingQTLsGTEx:
    message: 'Categorize splicing QTLs'
    output: 'data/CategorizeSplicingQTLs/GTEx/{tissue}.tsv'
    params:
        script = 'scripts/categorizeSqtls.R',
        qtl_files = 'code/results/qtl/noisy/GTEx/{tissue}/separateNoise/cis_100000/perm/chr*.addQval.txt.gz',
        fdr = 0.1,
    shell:
        '''
        Rscript {params.script} <(ls {params.qtl_files}) {params.fdr} {output}
        '''

rule GetEqtlDataGeuvadis:
    message: 'get eqtl data for Geuvadis' # This is from Ben. Notice there's a difference here, I selected 360 EUR individuals, but Ben's include both EUR and YRI (449 individuals)
    output:
        flag = touch('data/EQTLs/Geuvadis/done'),
        nominal = 'data/EQTLs/Geuvadis/steady_eqtl_nominal_frBen.txt.gz', # nominal file for coloc from Ben
        nominal_tbi = 'data/EQTLs/Geuvadis/steady_eqtl_nominal_frBen.txt.gz.tbi'
    params:
        nominal = '/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/QTLs/QTLTools/Expression.Splicing/NominalPassForColoc.txt.tabix.gz'
    shell:
        '''
        mkdir -p $(dirname {output.flag})
        cp {params.nominal} {output.nominal}
        cp {params.nominal}.tbi {output.nominal_tbi}

        echo "copied files to current directory: " > {output.flag}
        echo "cp {params.nominal} {output.nominal}" >> {output.flag}
        echo "cp {params.nominal}.tbi {output.nominal_tbi}" >> {output.flag}
        '''

use rule CategorizeSplicingQTLsGTEx as CategorizeSplicingQTLsGeuvadis with:
    output: 'data/CategorizeSplicingQTLs/Geuvadis/{population}.tsv'
    params:
        script = 'scripts/categorizeSqtls.R',
        qtl_files = 'code/results/qtl/noisy/Geuvadis/{population}/separateNoise/cis_100000/perm/chr*.addQval.txt.gz',
        fdr = 0.1
        
rule Ds_Dge_Results:
  output: 'data/ds_v_dge/{a_v_b}_data.rds'
  params:
    R_script = 'scripts/prepDGE_DS_AnalysesData.R',
    contrast = lambda w: w.a_v_b,
    gtf = '../annotations/hg38/gencode.v26.GRCh38.genes.csv',
    dsPrefix = 'code/results/ds/GTEx/',
    dgePrefix = 'code/results/dge/GTEx/',
    outPrefix = 'data/ds_v_dge/'
  shell:
    '''
    Rscript {params.R_script} {params.contrast} {params.gtf} {params.dsPrefix} {params.dgePrefix} {params.outPrefix}
    ls {output}
    '''





#WARN: Interactive rules!
#---------------------------------------------------------------
#WARN: Rerunning with snakeamke will overwrite ipynb!!!
# They are intended for creating jupyter notebooks, and using 
# the generated notebooks to interactively create plots and stuff.
#---------------------------------------------------------------                      

rule plot_DS_v_DGE:
#    run snakemake like so: 
#      snakemake --draft-notebook TARGET
#    then edit notebook in jupyter lab. Snakemake preemble will set up 
#    proper path and bring in input, output, params, etc.
  input: 
    rds_pat = expand('data/ds_v_dge/{avb}_data.rds', avb=CONTRASTS)
  output: 'plots/2024-03-08-plot_DS_v_DGE.pdf'
  notebook: 'scripts/2024-03-08-GTEx-DS-DGE-comps.r.ipynb'


