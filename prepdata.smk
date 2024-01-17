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





rule all:
    input: 
        # extract numerators from leafcutter2 with const run's results
        expand('data/ExtractFractions/GTEx/{tissue}.numerators_constcounts.noise_by_intron.rds',
               tissue=chosen_tissues),
        # categorize sQTLs into u-sQTLs and p-sQTLs, results in tsv file
        expand('data/CategorizeSplicingQTLs/GTEx/{tissue}.tsv',
               tissue=chosen_tissues),


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

use rule CategorizeSplicingQTLsGTEx as CategorizeSplicingQTLsGeuvadis with:
    output: 'data/CategorizeSplicingQTLs/Geuvadis/{population}.tsv'
    params:
        script = 'scripts/categorizeSqtls.R',
        qtl_files = 'code/results/qtl/noisy/Geuvadis/{population}/separateNoise/cis_100000/perm/chr*.addQval.txt.gz',
        fdr = 0.1,
        

