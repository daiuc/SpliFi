#---------------------------------------------------------------
#
# Author:      Chao Dai
# Description: Given summary stats, find GWAS hit loci given a pval
#              criteria, for instance 1e-7. See procedure illustration
#              below:
#              https://drive.google.com/file/d/1kAhHic0FqP2OPrl4STaELG2OdcugsC_v/view?usp=sharing
# Key inputs:  hg38 summary statistics of GWAS
# Note:        V2 does not reduce overlapping loci into a single merged locus.
#              Instead, it keeps overlapping locus as is.
#---------------------------------------------------------------




suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(GenomicRanges))
library(tictoc) # for fun


#-------------------------------------------------------------------------------
#                 SET UP RUN MODE
#-------------------------------------------------------------------------------

# RUNMODE:
# 1: snakemake
# 2: interactive or jupyter
# Note interactive() returns FALSE when run on jupyter
if (exists("snakemake")) RUNMODE = 1 else RUNMODE = 2

if (RUNMODE == 1) {
    sumstats.file = snakemake@input[['sumstats']]
    sumstats.colnames = snakemake@input[['colnames']]
    outHits.file = snakemake@output[[1]]
    pval.th = snakemake@params[['pval_th']]

    print("### Running in snakemake script mode.")

} else {
    sumstats.file = "resources/GWAS/Height_sumstats_hg38.tsv.gz"
    sumstats.colnames = "resources/GWAS/Height_sumstats_hg38.colnames.txt"
    outHits.file = NA
    pval.th = 1e-7
}

walk2(.x = c("input.summstats", "input.colnames", "output", "pval_threshold"),
      .y = c(sumstats.file, sumstats.colnames, outHits.file, pval.th),
      ~ cat(paste0("\n", .x, " : ", .y)))

#-------------------------------------------------------------------------------
#                 USEFUL FUNCTIONS
#-------------------------------------------------------------------------------


findGwasHitLocus = function(BagOfGranges, windowSize = 1e6, counter) {
#' @description : Find a locus, defined by `windowSize` flanking the topSNP
#'
#' @param BagofGranges : Remaining GenomicRanges (of snps) to be worked on
#' @param windowSize   : Window size of coloc, default 1e6
#' @param counter      :
#'
#' @return
#'      - HitLocus     : hit locus surrounding the lead SNP
#'      - RemainedSnps : remaining SNPs outside of the identified hit locus

    # find the leading SNP, and create 1mb window, center on SNP
    lead.snp = BagOfGranges[min(BagOfGranges$prank)]
    lead.snp = promoters(lead.snp,
                         upstream = as.integer(windowSize/2),
                         downstream = as.integer(windowSize/2))

    # remove all SNPs within locus from subsequent locus construct
    BagOfGranges = subsetByOverlaps(BagOfGranges, lead.snp,
                                    minoverlap = 1, invert = T)

    if (length(BagOfGranges) > 0) {
        # refresh pvalue rank
        BagOfGranges$prank = rank(BagOfGranges$prank, ties.method = "first")
    } else {
        print(paste0("Done. No more SNPs left after ", counter, " runs."))
    }

    return(list(HitLocus = lead.snp,
           RemainedSnps = BagOfGranges))
}


getNonOverlapLocus = function(x, locusSize = 1000000) {
  #' @description : compute non-overlapping windows.
  #' @param x     : a GRange object of snps that previously have overlapping
  #'                1mbp windows. x should all be on the same chromosome.
  #' @return      : a GRange object of redrawn locus windows.

  x = as.data.frame(sort(x)) %>%
    select(seqnames, pos = start, snp_id, pval) %>%
    mutate(d_prev = pos - lag(pos),
           d_next = lead(pos) - pos) %>%
    mutate_at(c("d_prev", "d_next"),
              ~ if_else(is.na(.x), as.integer(locusSize), .x)) %>%
    mutate(
      start = pos - as.integer(d_prev / 2) + 5,
      end = pos + as.integer(d_next / 2) - 5
    )
  x = makeGRangesFromDataFrame(x, keep.extra.columns = T, ignore.strand = T)
  return(x)
}



#-------------------------------------------------------------------------------
#               MAIN PROCEDURE
#-------------------------------------------------------------------------------


sumstats.colnames = read_lines(sumstats.colnames)
gwas = fread(sumstats.file, sep = "\t", col.names = sumstats.colnames)

# first remove any snps not passing threshold
gwas = gwas[pval < pval.th][order(pval)] %>% unique
gwas[, prank := rank(pval, ties.method = "first")]

print(paste0("Found ", nrow(gwas), " SNPs passing threshold."))



# Find Loci --------------------------------------------------------------------


# make GRanges
gwas.gr = makeGRangesFromDataFrame(gwas[, .(chr, pos, snp_id, pval, prank)],
                                   keep.extra.columns = T,
                                   ignore.strand = T,
                                   seqnames.field = "chr",
                                   start.field = "pos",
                                   end.field = "pos", )


# initialize
max.loop = length(gwas.gr)
counter = 1
locusSize = 1000000
loci = list()

# Loop over all the considered GWAS SNPs ton construct loci
while (length(gwas.gr) > 0 & counter < max.loop) {
    res = findGwasHitLocus(gwas.gr, windowSize = locusSize, counter = counter)
    loci = c(loci, res$HitLocus)
    gwas.gr = res$RemainedSnps
    counter = counter + 1
}

# convert list of granges to a single grange object and sort
loci = purrr::reduce(loci, c)
loci = sort(loci)
cat(paste0("Found ", length(loci), " GWAS hit loci.\n"))
cat("Loci may overlap.")



# Write output -----------------------------------------------------------------


as.data.table(loci) %>%
    .[,.(seqnames, start, end, snp_id, pval)] %>%
    write_tsv(outHits.file)

