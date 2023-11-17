

import gzip

# first parse annotations


transcripts = {}
for ln in gzip.open("gencode.v37.annotation.gff3.gz"):
    if ln[0] == "#": continue
    ln = ln.split()
    if ln[2] != "exon": continue
    chrom, s,e, strand, tname = ln[0], ln[3], ln[4], ln[6], ln[8].split("transcript_name=")[-1].split(";")[0]
    gname = ln[8].split("gene_name=")[-1].split(";")[0]

    tname = (gname,tname)
    if tname not in transcripts:
        transcripts[tname] = []

    s, e = int(s), int(e)
    transcripts[tname].append((chrom,s,e,strand))


introns = {}
intron2gene = {}
for t in transcripts:
    exons = transcripts[t]
    exons.sort()

    for i in range(len(exons)-1):
        #print t, exons[i],exons[i+1], exons[i][2],exons[i+1][1], exons[i][3]

        intron = (exons[i][0], int(exons[i][2]),int(exons[i+1][1]), exons[i][3])

        if intron not in introns:
            introns[intron] = 0
            intron2gene[intron] = {}
        introns[intron] += 1


        if t[0] not in intron2gene[intron]:
            intron2gene[intron][t[0]] = 0
        intron2gene[intron][t[0]] += 1


from misc_helper import stream_table
intron_used = {}
intron_class = {}
allintrons = set()

fout = gzip.open("intron_class.txt.gz",'w')

foutAnnot = gzip.open("intron_class_annot.txt.gz",'w')

N_func, N_noisy = 0, 0
for dic in stream_table(open("/project2/yangili1/yangili/snaptron_analysis/scripts/GTEX_usages_summary.csv"),'\t'):
    
    chrom,s,e,strand = dic['intron'].split(":")
    usage, N =  dic['atLeast5'].split('/')
   
    s, e = int(s)-1, int(e)
    if (chrom,s,e,strand) in intron2gene:
        annot = True
    else:
        annot = False

    intron_used[(chrom, s,e, strand)] = usage

    allintrons.add((chrom, s,e, strand))
    if annot:
        NF = max(intron2gene[(chrom,s,e,strand)].values())
    
    if annot:
        if NF >= 3:
            C = "putative_functional"
        else:
            if usage <= 0.1:
                C = "noisy"
            else:
                C = "putative_functional"
    else:
        if usage <= 0.1:
            C = "noisy"
        else:
            C = "putative_functional"

    fout.write("%s %d %d %s %s\n"%(chrom, int(s), int(e), strand, C))
    
    if annot:
        foutAnnot.write("%s %d %d %s %s\n"%(chrom, int(s), int(e), strand, "putative_functional"))
    else:
        foutAnnot.write("%s %d %d %s %s\n"%(chrom, int(s), int(e), strand, "noisy"))
    if C == "noisy":
        N_noisy += 1
    else:
        N_func += 1
    
    intron_class[(chrom,s,e,strand)] = C

    if annot: 
        intron2gene.pop((chrom,s,e,strand))

print "Number of functional: %d (%.3f%%) Number of noisy: %d"%(N_func, N_func/float(N_noisy)*100, N_noisy)

for intron in intron2gene:
    if max(intron2gene[intron].values()) >= 1:
        C = "putative_functional"
    else:
        C = "putative_noisy"

    chrom,s,e,strand = intron

    fout.write("%s %d %d %s %s\n"%(chrom, int(s), int(e), strand, C))


fout.close()
foutAnnot.close()



allintrons = list(allintrons.union(set(introns.keys())))
allintrons.sort()


fout = open("annot_by_usage.txt","w")

for intron in allintrons:
    try: p5 = str(intron_used[intron])
    except: p5 = "NA"
    
    try: c = intron_class[intron]
    except:
        c = "NA"

    if intron in introns:
        ann = introns[intron]
    else:
        ann = 0
        #print intron, 
    fout.write('%s %s %d\n'%(p5,c,  ann))

fout.close()
