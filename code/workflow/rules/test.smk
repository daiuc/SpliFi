# The goal of this script is to parse through a set of junction files to (i) cluster junctions, (ii) determine how many split-reads in that cluster are noisy according to functional intron annotation (intron_class.txt.gz).
# I would like you to test it, and also try to see whether you can come up with some sanity checks that it produces the right output (edited) 
# intron_class.txt.gz was generated using this script: /project2/yangili1/noisy_splicing_code/parse_noisy_metrics.py

# I would like you to follow the trail and see whether you can streamline the way I did this

# there are several steps and I could have made an error

# please let me know if you have any questions

rule test1:
    '''
    The goal of this script is to parse through a set of junction files to 
    (i) cluster junctions, 
    (ii) determine how many split-reads in that cluster are noisy according 
    to functional intron annotation (intron_class.txt.gz).
    
    I would like you to test it, and also try to see whether you can come up 
    with some sanity checks that it produces the right output (edited) 
    intron_class.txt.gz was generated using this script: 
    /project2/yangili1/noisy_splicing_code/parse_noisy_metrics.py
    
    I would like you to follow the trail and see whether you can streamline 
    the way I did this. there are several steps and I could have made an error

    '''
    input: 'resources/test.junclist.txt'
    output: touch('results/noisy/test.done') 
    params:
        out_prefix = 'results/noisy',
        intron_class = 'resources/intron_class.txt.gz'
    shell:
        '''
        python2 workflow/scripts/leafcutter_cluster_regtools_noisy_YIL.py \
            --juncfiles {input} \
            --outprefix {params.out_prefix} \
            --minclureads 10 \
            --minreads 3 \
            --mincluratio .001 \
            --includeconst `# include constitutive introns` \
            --noise {params.intron_class} \
            --rundir {params.out_prefix} \
            --offset 1
            #--cluster full_clusters_refined.reformat \
        '''