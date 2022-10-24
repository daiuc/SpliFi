# all misellaneous rules



# rule getconda:
#     '''get required conda environment. Only need to run once.'''
#     input: "workflow/envs/leafcutter.yaml"
#     output: touch("getconda.done")
#     conda:
#         "envs/leafcutter.yaml"
#     shell:
#         "echo done"


rule computeGTExJuncsHistogram:
    '''
    Compute histogram of counts for each junction among all 17000+ samples in GTEx.
    Processed data is used for some initial analysis as shown in 
        `analysis/validate-gtex-junc-file.pub.qmd`. 
    The point is to check if this downloaded GTEx junction file had filtered out
    low read counts.
    '''
    input: 'resources/GTEx/juncs/GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct.gz'
    output: 'resources/GTEx/juncs/GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions_histogram.txt'
    script: '../scripts/make_histogram_of_GTEx_juncs.py'