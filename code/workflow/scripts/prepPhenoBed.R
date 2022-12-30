#!/usr/bin/env Rscript


#' @description Convert a noisy splicing read count table to rank normalized following qtltools' bed file format
#' @author Chao Dai
#' 
#' 


suppressPackageStartupMessages(library("argparse"))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(RNOmni))
suppressMessages(library(GenomicRanges))



#---------------------------- Functions --------------------------------------


parseFracs = function(fracs) {
    #' parse fractions
    #' @param frac vectors of fraction strings eg. 1/2
    
    l = str_split(fracs, "/")
    num = map_int(l, ~ as.integer(.x[1]))
    denom = map_int(l, ~ as.integer(.x[2]))
    psi = map_dbl(l, ~ if_else(.x[2] == "0", 0, as.integer(.x[1])/as.integer(.x[2])))
    
    res = list(num = num, denom = denom, psi = psi)
    return(res)
}

getNSamplesPassMinNoisyReads = function(noisyreads, min_noisy_reads) {
    # return number of samples that pass noisy reads threshold
    # nreads: vector, number of noisy reads
    # threshold: number of samples having noisy reads > threshold
    
    check_threshold = noisyreads > min_noisy_reads
    num_samples_passed = sum(check_threshold)
    return(num_samples_passed)
}

filterNoisy = function(rawcounts, th_MeanClusterReads = 30,
                       th_MinSampleNumbers = 10,
                       th_MinNoisyReadsPerSample = 3) {

    dcols = names(rawcounts)[-1] # data cols
    
    # parsing PSI, num, and denom
    parsed = rawcounts[, -c("chrom")] %>% map(., parseFracs)
    
    # noisy reads
    nums = cbind(rawcounts[, .(chrom)], as.data.table(map(parsed, ~.x$num)))
    
    # cluster reads
   denoms = cbind(rawcounts[, .(chrom)],
        as.data.table(map(parsed, ~.x$denom)))
    
    # noisy psi
    psi = cbind(rawcounts[, .(chrom)],
        as.data.table(map(parsed, ~.x$psi)))

    # row index of noisy splicing that pass mean cluster reads filter
    idx_PassCluReads = denoms[, -c('chrom')] %>% 
        rowMeans  %>% `>`(th_MeanClusterReads) %>% which

    # row index of noisy splicing with at least x samples
    # having min y noisy reads
    idx_MinSamples = nums[, -c('chrom')] %>%
        apply(1, function(row) getNSamplesPassMinNoisyReads(
                    row, th_MinNoisyReadsPerSample)) %>%
        `>`(th_MinSampleNumbers)  %>%
        which

    psi = psi[intersect(idx_PassCluReads, idx_MinSamples)]
    nums = nums[intersect(idx_PassCluReads, idx_MinSamples)]
    denoms = denoms[intersect(idx_PassCluReads, idx_MinSamples)]
    return(list(psi=psi, nums = nums, denoms = denoms))
}

splitPhenoName = function(noisyname) {
    # noisyname : charactor vector, eg 'chr1:4692:5659:clu_3_-_*'
    # returns   : 
    splitnames = str_split(noisyname, "[:_]", simplify = T) %>% as.data.table
    splitnames = splitnames[, .(chrom = V1,
                    start = as.integer(V2),
                    end = as.integer(V3), 
                    strand = V6)]
    return(splitnames)
}

formatTable = function(dt, gtf) {
    # dt   : raw datatable, first column chrom: chr1:4692:5659:clu_3_-_*, subsequent table psi or counts per sample
    # gtf : GRange object with gene_name in mcols, annotation for gene name
    
    dcols = names(dt)[-1] # data cols
    splitname = splitPhenoName(dt$chrom)
    dt = dt[, c(splitname, list(pid = chrom), .SD), .SDcols = dcols]
    
    # annotate ranges by finding overlaps with anotation
    gr = makeGRangesFromDataFrame(dt[, 1:5], keep.extra.columns = T)
    hits = GenomicRanges::findOverlaps(gr, gtf, minoverlap = 50, ignore.strand = F)
    gr = gr[queryHits(hits)]
    mcols(gr)$gene_name = gtf[subjectHits(hits)]$gene_name
    gr = as.data.table(mcols(gr)) # select only pid and gene_name
    
    # only select 1 pid annotation in case of overlapping with 1+ genes
    gr = gr[, .(gene_name, R = rank(gene_name, ties.method = "first")), 
            by = "pid"][R == 1, .(pid, gid = gene_name)]
    
    # annotate gid for phenotypes that have a gene match
    dt = gr[dt, on = "pid"] # right outer join
    dt[is.na(gid), gid := pid] # replace gid=NA with gid=pid
    selcols = c("chrom", "start", "end", "pid", "gid", "strand", dcols)
    dt = dt[, ..selcols]
    
    return(dt)
    
}

normalizePSI = function(dt) {
    # ranknorm normalize psi
    # dt : filtered PSI data.table, following qtltools pheno bed format
    
    dcols = names(dt)[-c(1:6)] # data cols
    tcols = names(dt)[1:6] # text cols
    
    # first row normalization to N(0,1)
    dt.norm = dt[, ..dcols] %>% t %>% scale %>% t

    # then Inverse normal transform columns and round to 5 digits
    dt.norm = apply(dt.norm, 2, RankNorm) %>% apply(2, round, 5)
    
    # add back text cols
    dt.norm = cbind(dt[, ..tcols], dt.norm)
    
    return(dt.norm)
}


#---------------------------- Main --------------------------------------



parser <- ArgumentParser()

parser$add_argument("-I", "--inputfile", dest="inputfile",
        type="character", required=T,
        help="noisy splicing count file")

parser$add_argument("-A", "--anno", dest="anno",
        type="character", required=T,
        help="6 column bed format annotation file for gene names, no header.")

parser$add_argument("-O", "--outprefix", dest="outprefix",
        type="character", required=T,
        help="output qtltools compatible phenotype bed file")

parser$add_argument("-C", "--mincluster", dest="mincluster",
    type="integer", default=10,
    help="for each noisy splicing, must have on average C number of total reads for the cluster (default: 10)")

parser$add_argument("-S", "--minsample", dest="minsample",
    type="integer", default=10, 
    help="for each noisy splicing, at least S number of samples must have x (specified by -N) noisy splicing reads (default: 10)")

parser$add_argument("-N", "--minnoisy", dest="minnoisy",
    type="integer", default=3,
    help="for each noisy splicing, at least x (specified by -S) number of samples must have N noisy splicing reads (default: 3)")

args <- parser$parse_args()

input = args$inputfile
outprefix = args$outprefix
anno = args$anno
th_MeanClusterReads = args$mincluster # mean cluster reads threhold
th_MinSampleNumbers = args$minsample # number of samples having at least th_MinNoisyReadsPerSample noisy reads
th_MinNoisyReadsPerSample = args$minnoisy # at least th_MinSampleNumbers of samples having at least this amount of noisy reads


# read annotation file
# tryCatch({
# },  message = "Error! Make sure your annotation file is BED 6 column format.\n",
#     finally = q("no")
# )

cat(paste0("Read annotation: ", anno, "\n"))
gtf = fread(anno, sep = "\t", header = F,
        col.names = c("chrom", "start", "end", "gene_name", "score", "strand"))
gtf = makeGRangesFromDataFrame(gtf[, -c('score')], keep.extra.columns = T)

# read in noisy splicing count table
cat(paste0("Read raw count file: ", input, "\n"))
noisy = fread(input, sep = " ", header = TRUE)

th_MinSampleNumbers = min(as.integer((ncol(noisy)-1)/4), th_MinSampleNumbers)

# filter
cat(paste0(
    "Filter phenotypes with MinClusterReads=", th_MeanClusterReads,
    " and at least ", th_MinSampleNumbers, 
    " samples with min ", th_MinNoisyReadsPerSample, " reads.\n"
    ))
noisy = filterNoisy(noisy, th_MeanClusterReads, 
        th_MinSampleNumbers, th_MinNoisyReadsPerSample)


cat("Format table to be compatible with qtltools. \n")
# format psi table to be like qtltools phenotype bed
psi = formatTable(noisy$psi, gtf) 

cat("RankNorm normalize PSI. \n")
psi = normalizePSI(psi)

cols = names(psi)
cols = str_replace(cols, "chrom", "#Chr")
names(psi) = cols

chroms = unique(psi$`#Chr`)

for (ch in chroms) {
    outfile = paste0(outprefix, ".", ch, ".bed")
    dt = psi[`#Chr` == ch][order(start, end)]
    cat(paste0("Write phenotype to ", outfile, "\n"))
    fwrite(dt, outfile, sep = "\t")
}





