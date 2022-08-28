# https://github.com/griffithlab/regtools
# /home/yangili1/tools/regtools/build/regtools junctions extract -a 8 -i 50 -I 500000 bamfile.bam -o outfile.junc
# Using regtools speeds up the junction extraction step by an order of magnitude or more

import sys
import tempfile
import os
import gzip
import shutil


def main(options,libl):
    
    
    if options.cluster == None:
        pool_junc_reads(libl, options)
        refine_clusters(options)
        addlowusage(options)
    sort_junctions(libl, options)
    merge_junctions(options)
    get_numers(options)
    
    annotate_noisy(options)

def annotate_noisy(options):

    outPrefix = options.outprefix
    rundir = options.rundir    
    fnameout = "%s/%s"%(rundir,outPrefix)
    noisy_annotation = options.noiseclass

    dic_class = {"putative_functional":"F", "putative_noisy":"PN", "noisy":"N"}
    dic_usage = {}

    if noisy_annotation != None:

        if options.verbose:
            sys.stderr.write("Loading %s for noisy splicing classification..\n"%noisy_annotation)

        dic_noise = {}
        for ln in gzip.open(noisy_annotation):
            chrom,s,e,strand, classification = ln.split()
            dic_noise[(chrom,int(s),int(e)-1,strand)] = classification

        if options.verbose:
            sys.stderr.write("Loaded..\n")


    if not options.const:
        fname =  fnameout+"_perind.counts.gz"
    else:
        fname = fnameout+"_perind.constcounts.gz"

    noisyperind = fname.replace(".gz",".noise.gz")
    noisydiag = fname.replace(".gz",".noise_by_intron.gz")
    numers = fname.replace(".gz",".noise.gz").replace("perind",'perind_numers')

    fout = gzip.open(noisyperind,'w')
    foutdiag = gzip.open(noisydiag,'w')
    foutnumers = gzip.open(numers,'w')

    F = gzip.open(fname)
    ln = F.readline()
    fout.write(ln)
    foutdiag.write(ln)
    foutnumers.write(" ".join(ln.split()[1:])+'\n')

    noisy_clusters = {}
    clusters = {}

    current_cluster = None

    for ln in F:
        ln = ln.split()
        intron = ln[0]
        chrom, s, e, clu = intron.split(":")
        strand = clu.split("_")[-1]
    
        if current_cluster != clu:
            clusters[clu] = {"noise":[], "functional":[]}
            #print current_cluster, clu
            
            if current_cluster != None:
                #print [x[0] for x in clusters[current_cluster]['functional']], [x[0] for x in clusters[current_cluster]['noise']]
                if len(clusters[current_cluster]['noise']) > 0 and len(clusters[current_cluster]['functional']) > 0:
                    ID = "%s:%d:%d:%s"%(clusters[current_cluster]['functional'][0][0].split(":")[0], 
                                        min([int(x[0].split(":")[1]) for x in clusters[current_cluster]['functional']]), 
                                        max([int(x[0].split(":")[2]) for x in clusters[current_cluster]['functional']]), 
                                        clusters[current_cluster]['functional'][0][0].split(":")[3])

                    usages = [0 for i in range(len(clusters[current_cluster]['noise'][0])-1)]
                    totuse = [int(x.split("/")[1]) for x in clusters[current_cluster]['noise'][0][1:]]

                    for noise_intron in clusters[current_cluster]['noise']:
                        noise_use = [int(x.split("/")[0]) for x in noise_intron[1:]]
                        for i in range(len(noise_use)):
                            usages[i] += noise_use[i]
                    fout.write(ID+"_*" +" "+ " ".join(["%d/%d"%(usages[i],totuse[i]) for i in range(len(usages))])+'\n')
                    nID = tuple(ID.split(":"))
                    
                    foutnumers.write("%s:%s*:%s*:%s"%nID +" "+ " ".join(["%d"%(usages[i]) for i in range(len(usages))])+'\n')

                    #print ID+"_*" +" "+ " ".join(["%d/%d"%(usages[i],totuse[i]) for i in range(len(usages))])
                    for lnw in clusters[current_cluster]['functional']:
                        fout.write(" ".join(lnw)+'\n')
                        
                        foutnumers.write(lnw[0]+" "+ " ".join(["%s"%y.split('/')[0] for y in lnw[1:]])+'\n')
                
                clusters.pop(current_cluster)
            current_cluster = clu


        intronid = (chrom,int(s),int(e),strand)
        usages = [int(x.split("/")[0])/(float(x.split("/")[1])+0.1) for x in ln[1:]]
        if sum(usages) == 0: continue
        med = median(usages)


        if med >= 0.1:
            # If the median usage is above 0.1, consider it functional
            classification = "putative_functional"
        else:

            if intronid in dic_noise:
                # intron found in GTEx, so use its classification
                classification = dic_noise[intronid]
            else:
                # intron not found in GTEx or annotation and med < 0.1
                classification = "noisy"
            
        if classification not in dic_usage:
            dic_usage[classification] = []
        
        if dic_class[classification] in ["N","PN"]:
            noisy_clusters[clu] = ''
            clusters[clu]["noise"].append(ln)
        else:
            clusters[clu]["functional"].append(ln)

        foutdiag.write(intron+":%s"%dic_class[classification]+" "+' '.join(ln[1:])+'\n')


    if len(clusters[current_cluster]['noise']) > 0 and len(clusters[current_cluster]['functional']) > 0:
        ID = "%s:%d:%d:%s"%(clusters[current_cluster]['functional'][0][0].split(":")[0], 
                            min([int(x[0].split(":")[1]) for x in clusters[current_cluster]['functional']]), 
                            max([int(x[0].split(":")[2]) for x in clusters[current_cluster]['functional']]), 
                            clusters[current_cluster]['functional'][0][0].split(":")[3])
        
        usages = [0 for i in range(len(clusters[current_cluster]['noise'][0])-1)]
        totuse = [int(x.split("/")[1]) for x in clusters[current_cluster]['noise'][0][1:]]

        for noise_intron in clusters[current_cluster]['noise']:
            noise_use = [int(x.split("/")[0]) for x in noise_intron[1:]]
            for i in range(len(noise_use)):
                usages[i] += noise_use[i]
            fout.write(ID +" "+ " ".join(["%d/%d"%(usages[i],totuse[i]) for i in range(len(usages))])+'\n')

    foutdiag.close()
    fout.close()

    #for cl in dic_usage:                                                                                                                                                       
    #    print cl, median(dic_usage[cl]), len(dic_usage[cl])
    #print "%d/%d (%.2f%%) clusters with a noisy event..."%(len(noisy_clusters), len(clusters), float(len(noisy_clusters))/len(clusters)*100)

def addlowusage(options):
    global chromLst

    

    outPrefix = options.outprefix

    rundir = options.rundir
    pooled = "%s/%s_pooled"%(rundir,outPrefix)

    minclureads = int(options.minclureads)
    minreads = int(options.minreads)

    if options.cluster == None:
        refined_cluster = "%s/%s_refined"%(rundir,outPrefix)
    else:
        refined_cluster = options.cluster

    outFile = "%s/%s_refined_noisy"%(rundir,outPrefix)
    outFile_lowusageintrons = "%s/%s_lowusage_introns"%(rundir,outPrefix)
    fout = file(outFile,'w')
    fout_lowusage = file(outFile_lowusageintrons,'w')
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

def pool_junc_reads(flist, options):

    global chromLst

    outPrefix = options.outprefix
    rundir = options.rundir
    maxIntronLen = int(options.maxintronlen)
    checkchrom = options.checkchrom
    print("Max Intron Length:%d"%maxIntronLen)
    outFile = "%s/%s_pooled"%(rundir,outPrefix)
    
    by_chrom = {}
    for libl in flist:
        
        lib = libl.strip()
        if not os.path.isfile(lib):
            continue

        if options.verbose:
            sys.stderr.write("scanning %s...\n"%lib)

        if ".gz" in lib:
            F = gzip.open(lib)
        else:
            F = open(lib)

        for ln in F:
            
            lnsplit=ln.split()
            if len(lnsplit)<6: 
                sys.stderr.write("Error in %s \n" % lib)
                continue
        
            if len(lnsplit) == 12:
                chrom, A, B, dot, counts, strand, rA,rb, rgb, blockCount, blockSize, blockStarts = lnsplit
                if int(blockCount) > 2:  
                    print(ln, "ignored...")
                    continue
                Aoff, Boff = blockSize.split(",")[:2]
                A, B = int(A)+int(Aoff), int(B)-int(Boff)+1
            
            elif len(lnsplit) == 6:
                # old leafcutter junctions
                chrom, A, B, dot, counts, strand = lnsplit
                A, B = int(A), int(B)

            if checkchrom and (chrom not in chromLst): continue

            A, B = int(A), int(B)+int(options.offset)
    
            if B-A > int(maxIntronLen): continue                                                                                                                   
            try: by_chrom[(chrom,strand)][(A,B)] = int(counts) + by_chrom[(chrom,strand)][(A,B)]                                                                  
            except:                                                                                                                                                
                try: by_chrom[(chrom,strand)][(A,B)] = int(counts)                                                                                                 
                except: by_chrom[(chrom,strand)] = {(A,B):int(counts)}  


    fout = file(outFile, 'w')
    Ncluster = 0
    sys.stderr.write("Parsing...\n")
    for chrom in by_chrom:
        read_ks = [k for k,v in by_chrom[chrom].items() if v >= 3] # a junction must have at least 3 reads
	
	read_ks.sort()
        
        sys.stderr.write("%s:%s.."%chrom)
        if read_ks:
	    clu = cluster_intervals(read_ks)[0]
            for cl in clu:
    
                if len(cl) > 0: # 1 if cluster has more than one intron  
                    buf = '%s:%s '%chrom
                    for interval, count in [(x, by_chrom[chrom][x]) for x in cl]:
                        buf += "%d:%d" % interval + ":%d"%count+ " "
                    fout.write(buf+'\n')
                    Ncluster += 1
        sys.stderr.write("\nWrote %d clusters..\n"%Ncluster)
    fout.close()


def sort_junctions(libl, options):

    global chromLst


    noisy_annotation = options.noiseclass
    outPrefix = options.outprefix
    rundir = options.rundir
    checkchrom = options.checkchrom

    if options.cluster == None:
        refined_cluster = "%s/%s_refined_noisy"%(rundir,outPrefix)
        sys.stderr.write("Using %s as refined cluster...\n"%refined_cluster)
    else:
        refined_cluster = options.cluster

    runName = "%s/%s"%(rundir, outPrefix)

    exons, cluExons = {}, {}
    cluN = 0

    for ln in open(refined_cluster):
        chrom = ln.split()[0]
        cluN += 1
        for exon in ln.split()[1:]:
            A, B, count = exon.split(":")
            
            if chrom not in exons:
                exons[chrom] = {}
            if (int(A),int(B)) not in exons[chrom]:
                exons[chrom][(int(A),int(B))] = [cluN]
            else:
                exons[chrom][(int(A),int(B))].append(cluN)
            if cluN not in cluExons:
                cluExons[cluN] = []
            cluExons[cluN].append((chrom, A, B))

    merges = {}
    for ll in libl:
        lib=ll.rstrip()
        if not os.path.isfile(lib):
            continue
        libN = lib
        if libN not in merges:
            merges[libN] = []
        merges[libN].append(lib)

    fout_runlibs = file(runName+"_sortedlibs",'w')

    for libN in merges:
        libName = "%s/%s"%(rundir,libN.split('/')[-1])
        by_chrom = {}
        foutName = libName+'.%s.sorted.gz'%(runName.split("/")[-1])

        fout_runlibs.write(foutName+'\n')

        if options.verbose:   
            sys.stderr.write("Sorting %s..\n"%libN)
        if len(merges[libN]) > 1:
            if options.verbose:   
                sys.stderr.write("merging %s...\n"%(" ".join(merges[libN])))
        else:
            pass
        fout = gzip.open(foutName,'w')

        fout.write("chrom %s\n"%libN.split("/")[-1].split(".junc")[0])

        
        for lib in merges[libN]:
            if ".gz" in lib: F = gzip.open(lib)
            else: F = open(lib)

            for ln in F:

                lnsplit=ln.split()

                if len(lnsplit)<6:
                    sys.stderr.write("Error in %s \n" % lib)
                    continue

                if len(lnsplit) == 12:
                    chrom, A, B, dot, counts, strand, rA,rb, rgb, blockCount, blockSize, blockStarts = lnsplit
                    if int(blockCount) > 2:
                        print(ln, "ignored...")
                        continue
                    Aoff, Boff = blockSize.split(",")[:2]
                    A, B = int(A)+int(Aoff), int(B)-int(Boff)+1

                elif len(lnsplit) == 6:
                    # old leafcutter junctions                                                                                                                       
                    chrom, A, B, dot, counts, strand = lnsplit
                    A, B = int(A), int(B)
                    
                A, B = int(A), int(B) + int(options.offset)
            
                chrom = (chrom,strand)
                if chrom not in by_chrom:
                    by_chrom[chrom] = {}
                intron = (A, B)
                
                if intron in by_chrom[chrom]:
                    by_chrom[chrom][intron] += int(counts)
                else:
                    by_chrom[chrom][intron] = int(counts)
                
        for clu in cluExons:
            buf = []
            ks = cluExons[clu]
            ks.sort()
            #print clu, ks
            tot = 0
            usages = []
            for exon in ks:
                chrom, start, end = exon
                chrom = tuple(chrom.split(":"))
                start, end = int(start), int(end)

                if chrom not in by_chrom:
                    pass

                elif (start,end) in by_chrom[chrom]:
                    tot += by_chrom[chrom][(start,end)]

            for exon in ks:
                chrom, start, end = exon
                start, end = int(start), int(end)
                chrom = tuple(chrom.split(":"))
                chromID, strand = chrom

                intron = chromID, start, end+1, strand

                if chrom not in by_chrom:
                    buf.append("%s:%d:%d:clu_%d_%s 0/%d\n"%(chromID,start, end,clu, strand, tot))
                elif (start,end) in by_chrom[chrom]:
                    buf.append("%s:%d:%d:clu_%d_%s %d/%d\n"%(chromID,start, end, clu,strand, by_chrom[chrom][(start,end)], tot))
                else:
                    buf.append("%s:%d:%d:clu_%d_%s 0/%d\n"%(chromID,start, end,clu,strand, tot))
        
            fout.write("".join(buf))
        fout.close()
    fout_runlibs.close()

def refine_clusters(options):

    outPrefix = options.outprefix
    rundir = options.rundir
    minratio = float(options.mincluratio)
    minclureads = int(options.minclureads)
    minreads = int(options.minreads)

    inFile = "%s/%s_pooled"%(rundir,outPrefix)
    outFile = "%s/%s_refined"%(rundir,outPrefix)
    fout = file(outFile,'w')

    Ncl = 0
    for ln in open(inFile):
        clu = []
        totN = 0
        chrom = ln.split()[0]
        for ex in ln.split()[1:]:
            A, B, N = ex.split(":")
            clu.append(((int(A),int(B)), int(N)))
            totN += int(N)
        
        if totN < minclureads: continue

        if options.const:
            if len(clu) == 1:
                buf = '%s ' % chrom
                for interval, count in clu:
                    buf += "%d:%d" % interval + ":%d"%(count)+ " "
                Ncl += 1
                fout.write(buf+'\n')
        
        for cl in refine_linked(clu):            
            rc = refine_cluster(cl,minratio, minreads)
            if len(rc) > 0:
                for clu in rc:
                    buf = '%s ' % chrom
                    for interval, count in clu:
                        buf += "%d:%d" % interval + ":%d"%(count)+ " "
                    Ncl += 1
                    fout.write(buf+'\n')
    sys.stderr.write("Split into %s clusters...\n"%Ncl)
    fout.close()


def merge_junctions(options):    
    ''' function to merge junctions '''

    outPrefix = options.outprefix
    rundir = options.rundir
    
    fnameout = "%s/%s"%(rundir,outPrefix)

    flist = "%s/%s_sortedlibs"%(rundir, outPrefix)

    lsts = []
    for ln in open(flist):
        lsts.append(ln.strip())
    if options.verbose:
        sys.stderr.write("merging %d junction files...\n"%(len(lsts)))
    
    # Change 300 if max open file is < 300
    N = min([300, max([100, int(len(lsts)**(0.5))])])

    tmpfiles = []
    while len(lsts) > 1:    
        clst = []
        
        for i in range(0,(len(lsts)/N)+1): 
            lst = lsts[N*i:N*(i+1)]
            if len(lst) > 0:
                clst.append(lst)
        lsts = []
    
        for lst in clst:
            if len(lst) == 0: continue
            tmpfile = tempfile.mktemp()
            os.mkdir(tmpfile)
            foutname = tmpfile+"/tmpmerge.gz"
            fout = gzip.open(foutname,'w')
            
            merge_files(lst, fout, options)
            lsts.append(foutname)
            tmpfiles.append(foutname)
            fout.close()
            
    if not options.const:
        shutil.move(lsts[0], fnameout+"_perind.counts.gz")
    else:
        shutil.move(lsts[0], fnameout+"_perind.constcounts.gz")

def merge_files(fnames, fout, options):

    fopen = []
    for fname in fnames:
        if fname[-3:] == ".gz":
            fopen.append(gzip.open(fname))
        else:
            fopen.append(open(fname))

    finished = False
    N = 0
    while not finished:
        N += 1
        if N % 50000 == 0: 
            sys.stderr.write(".")
        buf = []
        for f in fopen:
            ln = f.readline().split()
            if len(ln) == 0: 
                finished = True
                break
            chrom = ln[0]
            data = ln[1:]
            if len(buf) == 0:
                buf.append(chrom)
            buf += data
        if len(buf) > 0:
            if buf[0] == "chrom":
                if options.verbose:
                    sys.stderr.write("merging %d files"%(len(buf)-1))
            fout.write(" ".join(buf)+'\n')
        else:
            break

    sys.stderr.write(" done.\n")
    for fin in fopen:
        fin.close()

def median(lst):
    n = len(lst)
    s = sorted(lst)
    return (sum(s[n//2-1:n//2+1])/2.0, s[n//2])[n % 2] if n else None


def cluster_intervals(E):
    ''' Clusters intervals together. '''
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
    '''
    Checks if A and B overlaps
    '''

    if A[1] < B[0] or B[1] < A[0]:
        return False
    else: return True

def refine_linked(clusters):

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


def refine_cluster(clu, cutoff, readcutoff):
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
        if (count/float(totN) >= cutoff and count >= readcutoff):
            intervals.append(inter)
            dic[inter] = count
        else:
            reCLU = True

    
    if len(intervals) == 0: return []
    
    # This makes sure that after trimming, the clusters are still good
    # after filtering.

    Atmp, B = cluster_intervals(intervals)
    A = []

    if len(Atmp) > 0:
        for cl in Atmp:
            if len(cl) == 1:
                A.append(cl)
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
        #else:
        #    return []
    NCs = []
    for c in A:
        if len(c) > 1:
            NC = refine_cluster([(x, dic[x]) for x in c], cutoff, readcutoff)
            NCs += NC
    #print "return", NCs
    return NCs


def get_numers(options):
    outPrefix = options.outprefix
    rundir = options.rundir

    if not options.const:                                                                                                                                                                                                                                                                                                                           
        fname = "%s/%s_perind.counts.gz"%(rundir,outPrefix)                                                                                                                                                                                                                                                                                      
        fnameout = "%s/%s_perind_numers.counts.gz"%(rundir,outPrefix)
    else:
        fname = "%s/%s_perind.constcounts.gz"%(rundir,outPrefix)
        fnameout = "%s/%s_perind_numers.constcounts.gz"%(rundir,outPrefix)
    
    input_file=gzip.open(fname, 'rb')
    fout = gzip.open(fnameout,'w')
    first_line=True
    
    for l in input_file:
        if first_line:
            fout.write(" ".join(l.strip().split(" ")[1:])+'\n') # print the sample names
            first_line=False
        else:
            l=l.strip()
            words=l.split(" ")            
            fout.write(words[0]+ " "+ " ".join( [ g.split("/")[0] for g in words[1:] ] ) +'\n')

    input_file.close()
    fout.close()

if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option("-j", "--juncfiles", dest="juncfiles",
                  help="text file with all junction files to be processed")

    parser.add_option("-o", "--outprefix", dest="outprefix", default = 'leafcutter',
                  help="output prefix (default leafcutter)")

    parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")

    parser.add_option("-r", "--rundir", dest="rundir", default='./',
                  help="write to directory (default ./)")
    
    parser.add_option("-l", "--maxintronlen", dest="maxintronlen", default = 100000,
                  help="maximum intron length in bp (default 100,000bp)")

    parser.add_option("-m", "--minclureads", dest="minclureads", default = 30,
                  help="minimum reads in a cluster (default 30 reads)")

    parser.add_option("-M", "--minreads", dest="minreads", default = 5,
                  help="minimum reads for a junction to be considered for clustering(default 5 reads)")

    parser.add_option("-p", "--mincluratio", dest="mincluratio", default = 0.001,
                  help="minimum fraction of reads in a cluster that support a junction (default 0.001)")

    parser.add_option("-c", "--cluster", dest="cluster", default = None,
                  help="refined cluster file when clusters are already made")

    parser.add_option("-k", "--checkchrom", dest="checkchrom", action="store_true",default = False,
                  help="check that the chromosomes are well formated e.g. chr1, chr2, ..., or 1, 2, ...")
    
    parser.add_option("-C", "--includeconst", dest="const", action="store_true",default = True,
                  help="also include constitutive introns")
    
    parser.add_option("-N", "--noise", dest="noiseclass", default = None,
                  help="Use provided intron_class.txt.gz to help identify noisy junction")

    parser.add_option("-f", "--offset", dest="offset", default = 0,
                  help="Offset sometimes useful for off by 1 annotations. (default 0)")

    (options, args) = parser.parse_args()

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
            sys.stderr.write("%s does not exist... check your junction files.\n"%junc)
            exit(0)
        libl.append(junc)

    chromLst = ["chr%d"%x for x in range(1,23)]+['chrX','chrY']+["%d"%x for x in range(1,23)]+['X','Y']

    main(options, libl)
