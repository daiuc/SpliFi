#!/usr/bin/env python

'''Identify noisy splicing

This script takes in junction files produced by regtools, then construct intron 
clusters. Next, constructed intron clusters are processed to identify rarely
spliced introns. Such rarely spliced introns are based on certain filtering
cut-offs. The result is a text file, following the same format as the standard 
leafcutter tool. The first column of each row indicates genome coordinates of
introns and the labeling of N[: noisy] or F[: Functional].

Junction files are processed using regtools:

    * https://github.com/griffithlab/regtools
    * /home/yangili1/tools/regtools/build/regtools junctions extract -a 8 -i \
      50 -I 500000 bamfile.bam -o outfile.junc
    * Using regtools speeds up the junction extraction step by an order of 
      magnitude or more

Procedure sequence: 

    1. pool_junc_reads
    2. refine_clusters
    3. addlowusage
    4. sort_junctions
    5. merge_junctions
    6. get_numers
    7. annotate_noisy

'''

import sys
import tempfile
import os
import gzip
import shutil


__author__    = "Yang Li, Chao Dai"
__email__     = "chaodai@uchicago.edu"
__status__    = "Development"
__version__   =  "0.1"


def cluster_intervals(E):
    '''Clusters intervals together
    
    Parameters:
    -----------
    E : list 
        list of tuples, e.g. [(start, end)]
    
    Returns: tuple
    ---------
    Eclusters : list of list
        Each element is a list of tuples (introns) clustered into 1 cluster
    cluster: list
        A list of tuples of introns
    '''
    E.sort()
    current = E[0]
    Eclusters, cluster = [], []

    i = 0
    while i < len(E):

        if overlaps(E[i], current):
            cluster.append(E[i])
        else:
            Eclusters.append(cluster)
            cluster = [E[i]]
        current = (E[i][0], max([current[1], E[i][1]]))
        i += 1

    if len(cluster) > 0:
        
        Eclusters.append(cluster)

    
    return Eclusters, E


def overlaps(A,B):
    '''Checks if A and B overlaps
    
    A : tuple
        start and end coordinates of 1 intron
    B : tuple
        start and end coordinates of another intron
    '''

    if A[1] < B[0] or B[1] < A[0]:
        return False
    else: return True


def pool_junc_reads(flist, options):
    '''Pool junction reads

    Parameters
    ----------
    flist : str
        The file list
    options : argparse object
        Passed arguments

    -------
    return
        No returns. Use side effects.

    ------
    Side Effects:
        write introns and counts by clusters. 
        format: [chrom]:[strand] [start]:[end]:[reads]
                e.g. chr17:+ 410646:413144:3 410646:413147:62
    '''

    global chromLst

    outPrefix = options.outprefix
    rundir = options.rundir
    maxIntronLen = int(options.maxintronlen)
    checkchrom = options.checkchrom
    print(f"Max Intron Length: {maxIntronLen}")
    outFile = f"{rundir}/{outPrefix}_pooled"
    
    # intron coordinates stored by chrom:strand
    # 2 level keys:
    #   level1 key: (chrom, strand)
    #   level2 key: (intron_start, intron_end)
    by_chrom = {}

    for libl in flist:
        lib = libl.strip()
        if not os.path.isfile(lib):
            continue

        if options.verbose:
            sys.stderr.write(f"scanning {lib}...\n")

        if ".gz" in lib:
            F = gzip.open(lib)
        else:
            F = open(lib)

        for ln in F:
            ln = ln.decode('utf-8') # convert bytes to string
            
            lnsplit=ln.split()
            if len(lnsplit) < 6: 
                sys.stderr.write(f"Error in {lib} \n")
                continue

            if len(lnsplit) == 12: # 12 fields regtools junc file
                chrom, A, B, dot, counts, strand, rA,rb, rgb, blockCount, blockSize, blockStarts = lnsplit
                if int(blockCount) > 2:  
                    print(ln, "ignored...")
                    continue
                Aoff, Boff = blockSize.split(",")[:2]
                A, B = int(A)+int(Aoff), int(B)-int(Boff)+1 # get intron

            elif len(lnsplit) == 6:
                # old leafcutter junctions
                chrom, A, B, dot, counts, strand = lnsplit
                A, B = int(A), int(B)
            
            if checkchrom and (chrom not in chromLst): 
                continue

            A, B = int(A), int(B)+int(options.offset)

            if B-A > int(maxIntronLen): continue
            
            # sum up all the reads at the same junctions if junctions already exist
            try: by_chrom[(chrom,strand)][(A,B)] = int(counts) + by_chrom[(chrom,strand)][(A,B)]                                                                  
            except:
                try: by_chrom[(chrom,strand)][(A,B)] = int(counts) # when only 1 junction
                except: by_chrom[(chrom,strand)] = {(A,B):int(counts)} # when only 1 junction


    with open(outFile, 'w') as fout:
        Ncluster = 0
        sys.stderr.write("Parsing...\n")
    
        for chrom in by_chrom:
            read_ks = [k for k,v in by_chrom[chrom].items() if v >= 3] # a junction must have at least 3 reads
            read_ks.sort()
        
            sys.stderr.write(f"{chrom[0]}:{chrom[1]}..")
            
            if read_ks:
                clu = cluster_intervals(read_ks)[0] # clusters of intervals
                
                for cl in clu:
                    if len(cl) > 0: # 1 if cluster has more than one intron  
                        buf = f'{chrom[0]}:{chrom[1]} '
                        for interval, count in [(x, by_chrom[chrom][x]) for x in cl]:
                            buf += f"{interval[0]}:{interval[1]}" + f":{count}" + " "
                        fout.write(buf+'\n')
                        Ncluster += 1
                        
        sys.stderr.write(f"\nWrote {Ncluster} clusters..\n")


def refine_linked(clusters):
    '''re-cluster introns into clusters of linked introns

    Linked introns are introns that share either 5' or 3' splice site
    
    Parameters:
    -----------
    clusters : tuple
        format is [((start, end), reads)], eg. 
        [((413430, 423479), 3), ((410646, 413144), 3), ((410646, 413147), 62), ((410646, 420671), 4)]
    
    Returns:
    --------
    return : list of list
        base element is a tuple, format: ((start, end), reads). e.g. 
        [[((413430, 423479), 3)], [((410646, 413144), 3), ((410646, 413147), 62), ((410646, 420671), 4)]]
    '''

    unassigned = [x for x in clusters[1:]]
    current = [clusters[0]]
    splicesites = set([current[0][0][0],current[0][0][1]]) # start & end of intron
    newClusters = []
    while len(unassigned) > 0:
        finished = False
    
        while not finished:
            finished = True
            torm = []
            for intron in unassigned:
                (start, end), count = intron
                if start in splicesites or end in splicesites:
                    current.append(intron)
                    splicesites.add(start)
                    splicesites.add(end)
                    finished = False
                    torm.append(intron)
            for intron in torm:
                unassigned.remove(intron)
        newClusters.append(current)
        current = []
        if len(unassigned) > 0:
            current = [unassigned[0]]
            splicesites = set([current[0][0][0],current[0][0][1]])
            unassigned = unassigned[1:]
    return newClusters


def refine_cluster(clu, cutoff, readcutoff):
    '''filter introns based on cutoffs

    Parameters:
    -----------
    clu : list of tuples, a single cluster
        list of tuples, each tuple an intron of the cluster
    cutoff : float
        reads ratio cutoff, passed in from option --mincluratio
    readcutoff : int
        minimum reads cutoff, passed in from option --minreads

    Filters:
    --------
        1. compute ratio of reads for each intron in a cluster
        2. remove intron if: 
            - ratio < ratio_cutoff
            - OR reads of the intron < readcutoff
        3. re-cluster with remaining introns
    
    Returns:
    --------
        return : list of list
        list of refined (filtered) clusters
    '''
    
    remove = []
    dic = {}
    intervals = []

    reCLU = False # re-cluster flag
    totN = 0

    for inter, count in clu:
        totN += count
    for inter, count in clu:
        if (count/float(totN) >= cutoff and count >= readcutoff):
            intervals.append(inter)
            dic[inter] = count # {(start, end): reads}
        else:
            reCLU = True # any intron not passing filters will enforce reCLU

    
    if len(intervals) == 0: return []
    
    # Below makes sure that after trimming/filtering, the clusters are still good
    # afterwards - each clusters have linked introns that pass filters.

    Atmp, _ = cluster_intervals(intervals)
    A = []

    # A: a list of linked intron clusters
    if len(Atmp) > 0:
        for cl in Atmp: # Atmp is a list of list
            if len(cl) == 1: # 1 intron
                A.append(cl)
            for c in refine_linked([(x,0) for x in cl]): # >1 introns
                if len(c) > 0:
                    A.append([x[0] for x in c])
            

    if len(A) == 1: # A has 1 cluster of introns
        rc = [(x, dic[x]) for x in A[0]]
    
        if len(rc) > 0:
            if reCLU: # recompute because ratio changed after removal of some introns
                return refine_cluster([(x, dic[x]) for x in A[0]], cutoff, readcutoff)
            else:
                return [[(x, dic[x]) for x in A[0]]]
        #else:
        #    return []
    
    NCs = [] # As in N Clusters, here A has more than 1 clusters of introns
    for c in A: 
        if len(c) > 1: # c has more than 1 introns
            NC = refine_cluster([(x, dic[x]) for x in c], cutoff, readcutoff)
            NCs += NC
    #print "return", NCs
    return NCs 



def refine_clusters(options):

    outPrefix = options.outprefix
    rundir = options.rundir
    minratio = float(options.mincluratio)
    minclureads = int(options.minclureads)
    minreads = int(options.minreads)

    inFile = f"{rundir}/{outPrefix}_pooled"
    outFile = f"{rundir}/{outPrefix}_refined"
    fout = open(outFile,'w')

    sys.stderr.write(f"Refine clusters from {inFile}...\n")

    Ncl = 0
    for ln in open(inFile): # pooled juncs
        clu = []
        totN = 0 # total cluster reads
        chrom = ln.split()[0]
        for ex in ln.split()[1:]:
            A, B, N = ex.split(":")
            clu.append(((int(A),int(B)), int(N)))
            totN += int(N)
        
        if totN < minclureads: continue

        if options.const: # include constitutive introns
            if len(clu) == 1:
                buf = f'{chrom} '
                for interval, count in clu:
                    buf += f"{interval[0]}:{interval[1]}" + f":{count}" + " "
                Ncl += 1
                fout.write(buf+'\n')
        
        for cl in refine_linked(clu):            
            rc = refine_cluster(cl,minratio, minreads)
            if len(rc) > 0:
                for clu in rc:
                    buf = f'{chrom} '
                    for interval, count in clu:
                        buf += f"{interval[0]}:{interval[1]}" + f":{count}" + " "
                    Ncl += 1
                    fout.write(buf+'\n')
    sys.stderr.write(f"Split into {Ncl} clusters...\n")
    fout.close()


# NOTE @ 22/9/5: still working on the addlowusage function
def addlowusage(options):
    '''
    
    Parameters:
    -----------
    options : argparse object
        pass in command options
    

    '''
    
    global chromLst

    
    outPrefix = options.outprefix
    rundir = options.rundir
    pooled = f"{rundir}/{outPrefix}_pooled"

    minclureads = int(options.minclureads)
    minreads = int(options.minreads)

    if options.cluster == None:
        refined_cluster = f"{rundir}/{outPrefix}_refined"
    else:
        refined_cluster = options.cluster

    outFile = f"{rundir}/{outPrefix}_refined_noisy"
    outFile_lowusageintrons = f"{rundir}/{outPrefix}_lowusage_introns"

    fout = open(outFile,'w')
    fout_lowusage = open(outFile_lowusageintrons,'w')
    
    exons5,exons3, cluExons = {}, {}, {}
    cluN = 0

    for ln in open(refined_cluster):
        chrom = ln.split()[0]
        cluN += 1
        cluExons[(chrom,cluN)] = []
        for exon in ln.split()[1:]:
            A, B, count = exon.split(":")
            if chrom not in exons5:
                exons5[chrom] = {}
                exons3[chrom] = {}
            exons5[chrom][int(A)] = (chrom,cluN)
            exons3[chrom][int(B)] = (chrom,cluN)
            cluExons[(chrom,cluN)].append(exon)
    
    lowusage_intron = {}
    test = ""
    
    for ln in open(pooled):
        #print "!", ln
        clu = []
        totN = 0
        chrom = ln.split()[0]
	#print chrom
        if chrom in exons5:
        
	    for exon in ln.split()[1:]:
                A, B, N = exon.split(":")
                
                if int(A) in exons5[chrom]:
                    clu = exons5[chrom][int(A)]
                    if exon not in cluExons[clu]:
                        cluExons[clu].append(exon)
                        if clu not in lowusage_intron:
                            lowusage_intron[clu] = []
                        lowusage_intron[clu].append(exon)

                elif int(B) in exons3[chrom]:
                    clu = exons3[chrom][int(B)]
                    if exon not in cluExons[clu]:
                        cluExons[clu].append(exon)
                        if clu not in lowusage_intron:
                            lowusage_intron[clu] = []
                        lowusage_intron[clu].append(exon)

                else:# int(A) not in exons5[chrom] and int(B) not in exons3[chrom]:
                    if int(N) > minreads:
                        cluN += 1
                        cluExons[(chrom, cluN)] = [exon]

    ks = lowusage_intron.keys()
    ks.sort()
    
    for clu in ks:
        fout_lowusage.write(clu[0] + " " + " ".join(lowusage_intron[clu])+'\n')
    fout_lowusage.close()


    cluLst = cluExons.keys()
    cluLst.sort()
    
    for clu in cluLst:
        if not options.const:
            if len(cluExons[clu]) == 1: continue

        if sum([int(ex.split(":")[-1]) for ex in cluExons[clu]]) < minclureads:
            continue
        chrom = clu[0]
        buf = '%s' % chrom
        for ex in cluExons[clu]:
            buf += " " + ex
        fout.write(buf+'\n')
    fout.close()



def main(options, libl):
    if options.cluster == None:
        pool_junc_reads(libl, options)
        refine_clusters(options)


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-j", "--juncfiles", dest="juncfiles",
        type=str, required=True,
        help="text file with all junction files to be processed")

    parser.add_argument("-o", "--outprefix", dest="outprefix",
        default = 'leafcutter',
        help="output prefix (default leafcutter)")

    parser.add_argument("-q", "--quiet", dest="verbose", default=True,
        action="store_false", help="don't print status messages to stdout")

    parser.add_argument("-r", "--rundir", dest="rundir", default='./',
                  help="write to directory (default ./)")
    
    parser.add_argument("-l", "--maxintronlen", dest="maxintronlen",
        default = 100000, 
        help="maximum intron length in bp (default 100,000bp)")

    parser.add_argument("-m", "--minclureads", dest="minclureads", default = 30,
        help="minimum reads in a cluster (default 30 reads)")

    parser.add_argument("-M", "--minreads", dest="minreads", default = 5,
                  help="minimum reads for a junction to be considered for \
                        clustering(default 5 reads)")

    parser.add_argument("-p", "--mincluratio", dest="mincluratio", 
        default = 0.001,
        help="minimum fraction of reads in a cluster that support a junction \
              (default 0.001)")

    parser.add_argument("-c", "--cluster", dest="cluster", default = None,
        help="refined cluster file when clusters are already made")

    parser.add_argument("-k", "--checkchrom", dest="checkchrom",
        action="store_true",default = False,
        help="check that the chromosomes are well formated e.g. chr1, chr2, \
              ..., or 1, 2, ...")
    
    parser.add_argument("-C", "--includeconst", dest="const", \
        action="store_true", default = True, 
        help="also include constitutive introns")
    
    parser.add_argument("-N", "--noise", dest="noiseclass", default = None,
        help="Use provided intron_class.txt.gz to help identify noisy junction")

    parser.add_argument("-f", "--offset", dest="offset", default = 0,
        help="Offset sometimes useful for off by 1 annotations. (default 0)")

    options = parser.parse_args()

    if options.juncfiles == None:
        sys.stderr.write("Error: no junction file provided...\n")
        exit(0)
    
    # Get the junction file list
    libl = []
    for junc in open(options.juncfiles):
        junc = junc.strip()
        try:
            open(junc)
        except: 
            sys.stderr.write(f"{junc} does not exist... check your junction files.\n")
            exit(0)
        libl.append(junc)

    chromLst = [f"chr{x}" for x in range(1,23)]+['chrX','chrY']+[f"{x}" for x in range(1,23)]+['X','Y']

    main(options, libl)
