
import sys
import os
import gzip
# import shutil


def aggJuncs(juncfile, clusterfile, df, rundir, offset):
    '''Aggregate junction files by group / individual

    In Geuvadis data, an individual ID may cooresponds to multiple run IDs, eg.
    1 NA##### -> multiple ERR###### files. This function aggregates these 
    individual ERR* files into a single per individual NA* file. 

    Parameters:
    -----------
        juncfile    : a text file
            A list of junction file paths
        clusterfile : a text file
            A file for clustered introns
        df      :  dataframe
            Dataframe with SUBJID (subject ID), eg. the NA#### in Geuvadis
            and SAMPID (sample ID), eg. the ERR###### in Geuvadis. Both 
            columns are also in index.
        rundir  : str
            run directory
        offset  : int
    
    Returns:
    --------
        return : no returns. Use site effect.

    Side-effects:
    -------------
        text file:  '{rundir}/subjid.agg.junc.gz'
        merged and sorted junction file with read counts aggregated per individual. 
    '''

    with open(juncfile) as f:
        libl = [ln.strip() for ln in f.readlines()] # list of junc file paths
    

    #--------- clusters (provided) -------------
    with open(clusterfile) as f:
        # exons:  { k=chrom : v={ k=(start, end) : v=clusterID } }
        # cluExons:  { k=clusterID : v=[(chrom, start, end)] }
        exons, cluExons = {}, {} 
        cluN = 0 # clusterID
        
        # fill in exons, cluExons dict from `*refined_noisy` intron cluster file
        for ln in f: # "chr10:+ 135203:179993:5 135302:179993:29"
            chrom = ln.split()[0] # e.g. "chr10:+"
            cluN += 1
            for exon in ln.split()[1:]: # e.g. "135203:179993:5 135302:179993:29"
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

    #--------- unmerged junc files -------------
    merges = {} # stores junc file names as dict { k=subject : v=[filenames] }
    for lib in libl:
        # here libN = SUBJID; lib is same level as sample ID
        sampN = lib.split('/')[-1].split('.')[0] # sample name eg. GTEX-1117F-0226-SM-5GZZ7
        subjN = df.query('SAMPID == @sampN')['SUBJID'][0] # subject id, eg. GTEX-1117F

        if not os.path.isfile(lib):
            continue
        if subjN not in merges:
            merges[subjN] = []
        merges[subjN].append(lib)

    # operate on each subject, each subjN can have more than 1+ junc files
    for subjN in merges: 
        by_chrom = {} # to store junctions from original unsorted junc file
        
        
        # write sorted junc file names into intermediate file
        foutName = os.path.join(rundir, subjN + '.agg.junc.gz') # 'test/gtex_w_clu/GTEX-1IDJU.agg.junc.gz'

        sys.stderr.write(f"Merging {subjN}..\n")
        if len(merges[subjN]) > 1: 
            sys.stderr.write(f"merging {' '.join(merges[subjN])}...\n")
        else: pass
        fout = gzip.open(foutName,'wt') 

        #-------- Process and write junction files --------

        #-------- Gather counts from all junc files of library --------
        # store in by_chrom: { ('chr1', '+') : { (100, 300) : 5, (500, 700): 10, ... } }
        for lib in merges[subjN]:
            if ".gz" in lib: 
                F = gzip.open(lib)
            else: 
                F = open(lib)
        
            for ln in F: # 1 line: e.g. "chr17\t81701131\t81701534\t.\t1\t+"

                if type(ln) == bytes:
                    ln = ln.decode('utf-8') # convert bytes to string
                lnsplit = ln.split()

                if len(lnsplit) < 6:
                    sys.stderr.write(f"Error in {lib} \n")
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
                    
                A, B = A, B + int(offset) # start, end + offset
            
                chrom = (chrom, strand)
                if chrom not in by_chrom: 
                    by_chrom[chrom] = {} # store introns from junc file, key: ('chr1', '+')
                
                intron = (A, B)
                if intron in by_chrom[chrom]: # sum up reads by intron from junc files
                    by_chrom[chrom][intron] += int(counts)
                else:
                    by_chrom[chrom][intron] = int(counts)
        
        
        #------- Take clusters from refined_noisy, assign reads -------
        # reads are from by_chrom (junc files)
        # For each intron cluster, write fraction for each intron (one intron per line).
        for clu in cluExons: # cluExons: { k=cluID : v=[(chrom, start, end)...]}
            buf = []
            ks = cluExons[clu] # eg: [('chr1:+', 827776, 829002), ..] introns of a clu
            ks.sort() # no need to version sort within cluster

        #     # Step 1: sum cluster level reads from each intron
        #     # gather (sum) reads for each cluster in refined_noisy, read counts are from junc file (by_chrom)
        #     tot = 0 # sum of total read counts per cluster
        #     usages = []
        #     for exon in ks:
        #         chrom, start, end = exon
        #         chrom = tuple(chrom.split(":")) # convert 'chr3:+' to ('chr3', '+') as in by_chrom
        #         start, end = int(start), int(end)

        #         if chrom not in by_chrom:
        #             pass
        #         elif (start, end) in by_chrom[chrom]:
        #             tot += by_chrom[chrom][(start,end)]

            # Step 2: append intron usage fraction to stream buffer
            for exon in ks:
                chrom, start, end = exon
                start, end = int(start), int(end)
                chrom = tuple(chrom.split(":"))
                chromID, strand = chrom # chromID eg: 'chr3'

                intron = chromID, start, end+1, strand # converting to 1-based coordinates

                if chrom not in by_chrom: 
                    # if refined exon chrom is not found in junc file, write 0
                    #buf.append(f"{chromID}\t{start}\t{end}\t.\t0\t{strand}\n")
                    pass

                elif (start,end) in by_chrom[chrom]: 
                    # if refind exon is in junc file, write exon reads 
                    buf.append(f"{chromID}\t{start}\t{end}\t.\t{by_chrom[chrom][(start,end)]}\t{strand}\n")
                else:
                    # if refined exon is not found in junc file, write 0
                    # buf.append(f"{chromID}\t{start}\t{end}\t.\t0\t{strand}\n")
                    pass
                
        
            fout.write("".join(buf))
        fout.close()


def main(options, df):
    aggJuncs(juncfile=options.juncfiles, clusterfile=options.cluster,
             df=df, rundir=options.rundir, offset=options.offset)

#------------------------------------------------------------------

if __name__ == '__main__':

    import argparse
    import pandas as pd

    parser = argparse.ArgumentParser()

    parser.add_argument("-j", "--juncfiles", dest="juncfiles",
        type=str, required=True,
        help="text file with all junction files to be processed")

    parser.add_argument("-c", "--cluster", dest="cluster", default = None,
        help="refined cluster file when clusters are already made")
    
    parser.add_argument("-L", "--lookup", dest="lookup", default = None,
        help="lookup tsv, first col SAMPID, second col SUBJID")

    parser.add_argument("-f", "--offset", dest="offset", default = 0,
        help="Offset sometimes useful for off by 1 annotations. (default 0)")

    parser.add_argument("-r", "--rundir", dest="rundir", default='./',
                help="write to directory (default ./)")
    
    options = parser.parse_args()

    if options.juncfiles == None:
        sys.stderr.write("Error: no junction file provided...\n")
        exit(0)
    
    if options.lookup == None:
        sys.stderr.write("Error: no lookup table provided...\n")
        exit(0)
    else:
        df = pd.read_csv(options.lookup, sep='\t')
        if list(df.columns) == ['SAMPID', 'SUBJID']:
            df.set_index(['SAMPID', 'SUBJID'], drop=False, inplace=True)
        else:
            sys.stderr.write("Error: Dataframe must be strictly col1=SAMPID, col2=SUBJID...\n")
            exit(0)
    
    main(options, df)