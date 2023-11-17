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
suppressMessages(library(furrr))




#---------------------------- Functions --------------------------------------


# improved parse function
parseFracs = function(fracs) {
  #' parse fractions
  #' @param frac vectors of fraction strings eg. c("1/2", "2/3")

  df <- data.frame(f = fracs)
  df <- separate(df, col = f, into = c("num", "denom"), sep = "/", convert = T)
  df <- mutate(df, psi = num / (denom + 1e-3))
  return(df)
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

    keep_row_ids <- intersect(idx_PassCluReads, idx_MinSamples)
    psi = psi[keep_row_ids]
    nums = nums[keep_row_ids]
    denoms = denoms[keep_row_ids]
    return(list(psi=psi, nums = nums, denoms = denoms))
}


splitPhenoName = function(noisyname) {
  # noisyname : charactor vector, eg 'chr1:4692:5659:clu_3_-_*'
  # returns   : dataframe of chrom, start, end, strand columns

  df <- data.frame(v = noisyname)
  df <- extract(df, v, into = c("chrom", "start", "end", "strand"),
                regex = "(chr[0-9XY]{1,2}):(\\d+):(\\d+):clu[_\\d]+_([\\+\\-])",
                convert = T)

  return(df)
}



formatTable = function(dt, gtf) {
    # dt   : raw datatable, first column chrom: chr1:4692:5659:clu_3_-_*, subsequent table psi or counts per sample
    # gtf : GRange object with gene_id in mcols

    dcols = names(dt)[-1] # data cols
    splitname = splitPhenoName(dt$chrom)
    dt = dt[, c(splitname, list(pid = chrom), .SD), .SDcols = dcols]

    # annotate ranges by finding overlaps with anotation
    gr = makeGRangesFromDataFrame(dt[, 1:5], keep.extra.columns = T)
    # hits = GenomicRanges::findOverlaps(gr, gtf, minoverlap = 50, ignore.strand = F)
    hits = GenomicRanges::findOverlaps(gr, gtf, minoverlap = 10, ignore.strand = F)
    gr = gr[queryHits(hits)]
    mcols(gr)$gene_id = gtf[subjectHits(hits)]$gene_id
    gr = as.data.table(mcols(gr)) # select only pid and gene_id

    # only select 1 pid annotation in case of overlapping with 1+ genes
    gr = gr[, .(gene_id, R = rank(gene_id, ties.method = "first")),
            by = "pid"][R == 1, .(pid, gid = gene_id)]

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
    dt.norm = dt[, ..dcols] %>% t() %>% scale() %>% t()

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

parser$add_argument("-T", "--threads", dest="threads",
    type="integer", default=1,
    help="number of threads")

args <- parser$parse_args()

# multi session
nCores <- min(availableCores(), args$threads)
plan(strategy = multisession, workers = nCores)

# other inputs
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
gtf = fread(anno) %>%
  .[feature == "gene" & gene_type == "protein_coding",
    .(seqname, start, end, strand, gene_id)]
gtf <- unique(gtf)
gtf = makeGRangesFromDataFrame(gtf, keep.extra.columns = T)

# read in noisy splicing count table
cat(paste0("Read raw count file: ", input, "\n"))
noisy = fread(input, sep = " ", header = TRUE)

# determine min sample numbers satisfying min reads
th_MinSampleNumbers = min(as.integer((ncol(noisy)-1)/4), th_MinSampleNumbers)

# split count table by chromosome

# Create a grouping variable
noisy <- extract(noisy, chrom, "group", "(chr[\\dXY]{1,2})", remove = F) %>% data.table()

# Split the data.table based on the grouping variable
noisy_split <- split(noisy, by = "group", keep.by = FALSE)

# remove group column from original dt
noisy[, group := NULL]




# filter
cat(paste0(
    "Filter phenotypes with MinClusterReads=", th_MeanClusterReads,
    " and at least ", th_MinSampleNumbers,
    " samples with min ", th_MinNoisyReadsPerSample, " reads.\n"
    ))

noisy_split <- future_map(noisy_split,
  ~filterNoisy(.x,
               th_MeanClusterReads,
               th_MinSampleNumbers,
               th_MinNoisyReadsPerSample))


cat("Format table to be compatible with qtltools. \n")

# format psi table to be like qtltools phenotype bed
psi_split = future_map(noisy_split, ~formatTable(.x$psi, gtf))



cat("RankNorm normalize PSI genomewide. \n")
# psi_split = future_map(psi_split, ~normalizePSI(.x))
psi_dt = rbindlist(psi_split)
psi_dt = normalizePSI(psi_dt)

# split by chrom
psi_split = split(psi_dt, by = "chrom", keep.by = TRUE)

# sort list names
psi_split <- psi_split[gtools::mixedsort(names(psi_split))]


# write results
cols = names(psi_split[[1]]) %>% str_replace(., "chrom", "#Chr")

for (ch in names(psi_split)) {
  outfile = paste0(outprefix, ".", ch, ".bed")
  dt = psi_split[[ch]]
  names(dt) = cols
  dt = dt[`#Chr` == ch][order(start, end)]
  cat(paste0("\nWrite phenotype to ", outfile, "\n"))
  fwrite(dt, outfile, sep = "\t")
  cat(glue::glue("\n{outfile} written.\n"))
}




#-------------------- test: comment out ----------------------------------------
#
# suppressMessages(library(furrr))
#
#
# nCores <- availableCores()
# plan(strategy = multisession, workers = nCores)



# anno <- "/project2/yangili1/cdai/genome_index/hs38/gencode.v38.primary_assembly.annotation.dataframe.csv"
# input <- "results/pheno/noisy/Geuvadis/EUR/leafcutter_perind.counts.inclProductive.gz"
# th_MeanClusterReads <- 10
# th_MinSampleNumbers <- 3
# th_MinNoisyReadsPerSample <- 2
#
#
# noisy <- noisy[sample(1:nrow(noisy), 5000)]




