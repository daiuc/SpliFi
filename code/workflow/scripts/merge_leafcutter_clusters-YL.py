import sys

try: 
    name = sys.argv[1]
    minread = sys.argv[2]
    clusters = sys.argv[3:]
except:
    sys.stderr.write("Usage: python merge_leafcutter_cluster.py combined_refined 20 refined_clusters_1 ...\n")
    exit()

try: open(name)
except:
    pass
else:
    sys.stderr.write("%s exists already..\nrm %s to continue"%(name,name))
    #exit()
fout = file(sys.argv[1],'w')
    
def main(minread):
    by_chrom = {}
    N = 0
    for cl in clusters:
        print "%s %d of %d (%.2f%%)"%(cl, N, len(clusters), 100*float(N)/len(clusters))
        N += 1
        for ln in open(cl):
            ln = ln.split()
            chrom = tuple(ln[0].split(":"))
            chrom = ("chr"+chrom[0].replace("chr",""), chrom[1])
            if chrom not in by_chrom:
                by_chrom[chrom] = {}
            junctions = ln[1:]
            for j in junctions:
                s,e, c = j.split(':')
                if (int(s),int(e)) not in by_chrom[chrom]:
                    by_chrom[chrom][(int(s),int(e))] = 0
                by_chrom[chrom][(int(s),int(e))]+= int(c)

    for (chrom,strand) in by_chrom:
        read_ks = [k for k,v in by_chrom[(chrom,strand)].items()]
        read_ks.sort()
        clu = cluster_intervals(read_ks)[0]
        for cl in clu:
            #if len(cl) > 1: # if cluster has more than one intron
        
            buf = '%s:%s '%(chrom,strand)
            cluster = [(x, by_chrom[(chrom,strand)][x]) for x in cl]
            #print chrom, strand, cluster
            #print "---\n", cl
            
            for fclu in refine_clusters(cluster, minread):
                fout.write(" ".join(["%s:%s"%(chrom,strand)]+ ["%s:%d"%("%s:%s"%x[0],x[1]) for x in fclu])+'\n')
 
    fout.close()

def refine_clusters(clu, minread):
    if True:
        Ncl = 0
        for cl in refine_linked(clu):
            rc = refine_cluster(cl, readcutoff=minread)
            if len(rc) > 0:
                for clu in rc:
                    #buf = '%s ' % chrom
                    #for interval, count in clu:
                    #    buf += "%d:%d" % interval + ":%d"%(count)+ " "
                    Ncl += 1
                    yield clu
    #sys.stderr.write("Split into %s clusters...\n"%Ncl)
    

def cluster_intervals(E):
    ''' Clusters intervals together. '''
    E.sort()
    if len(E) == 0:
        return [], []
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
    '''
    Checks if A and B overlaps
    '''

    if A[1] < B[0] or B[1] < A[0]:
        return False
    else: return True



def refine_linked(clusters):
    if len(clusters) == 1:
        return [clusters]
    unassigned = [x for x in clusters[1:]]
    current = [clusters[0]]
    splicesites = set([current[0][0][0],current[0][0][1]])
    newClusters = []
    while len(unassigned) > 0:
        finished = False
    
        while not finished:
            finished = True
            torm = []
            for intron in unassigned:
                inter, count = intron
                start, end = inter
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

def refine_cluster(clu, readcutoff = 20, cutoff=0.01):
    ''' for each exon in the cluster compute the ratio of reads, if smaller than cutoff,
    remove and recluster '''
    
    remove = []
    dic = {}
    intervals = []

    reCLU = False
    totN = 0

    for inter, count in clu:
        totN += count
    for inter, count in clu:
        if ((count/float(totN)) >= cutoff and (count >= 10)) or (count >= int(readcutoff)):
            intervals.append(inter)
            dic[inter] = count
        else:
            reCLU = True
            
    if len(intervals) == 0: return []
    Atmp, B = cluster_intervals(intervals)
    A = []
    for cl in Atmp:
        for c in refine_linked([(x,0) for x in cl]):
            if len(c) > 0:
                A.append([x[0] for x in c])
    
    if len(A) == 1:
        rc = [(x, dic[x]) for x in A[0]]
        if len(rc) > 0:
            if reCLU:
                return refine_cluster([(x, dic[x]) for x in A[0]], cutoff, readcutoff)
            else:
                return [[(x, dic[x]) for x in A[0]]]
        else:
            return []
    NCs = []
    for c in A:
        if len(c) > 0:
            NC = refine_cluster([(x, dic[x]) for x in c], cutoff, readcutoff)
            NCs += NC
    return NCs



main(minread)
