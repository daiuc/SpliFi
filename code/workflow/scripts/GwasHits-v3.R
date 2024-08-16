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
# 
#              This is copied and modifed from the A2I project script.
#---------------------------------------------------------------



if (interactive()) {
    args = scan(
                text = "
                resources/GWAS/hg38/HT.tsv.gz 
                1e-5 
                ", what = character())
    sumstats.file = args[1]
    pval.th = as.numeric(args[2])
} else {
    args = commandArgs(trailingOnly = TRUE)
    if (length(args) != 3) {
        stop("Usage: Rscript script.R sumstats.file outHits.file pval_threshold")
    }
    sumstats.file = args[1]
    outHits.file = args[2]
    pval.th = as.numeric(args[3])
}



suppressMessages(library(tidyverse))
suppressMessages(library(glue))
suppressMessages(library(data.table))
suppressMessages(library(GenomicRanges))



#-------------------------------------------------------------------------------
#                  FUNCTIONS
#-------------------------------------------------------------------------------


findGwasHitLocus = function(BagOfGranges, locusWidth = 1e6) {
#' @description : Find a locus, defined by `locusWidth` flanking the topSNP
#'
#' @param BagofGranges : Remaining GenomicRanges (of snps) to be worked on
#' @param locusWidth   : Window size of coloc, default 1e6
#'
#' @return
#'      - HitLocus     : hit locus surrounding the lead SNP
#'      - RemainedSnps : remaining SNPs outside of the identified hit locus

    # find the leading SNP, and create 1mb window, center on SNP
    minSnpRank = min(BagOfGranges$prank)
    lead.snp = BagOfGranges[BagOfGranges$prank == minSnpRank]
    lead.snp = promoters(lead.snp,
                 upstream = as.integer(locusWidth/2),
                 downstream = as.integer(locusWidth/2))

    # remove all SNPs within locus from subsequent locus construct
    BagOfGranges = subsetByOverlaps(
                        BagOfGranges,
                        lead.snp,
                        minoverlap = 1,
                        invert = T)

    return(list(HitLocus = lead.snp, RemainedSnps = BagOfGranges))
}


#-------------------------------------------------------------------------------
#               LOAD FILES
#-------------------------------------------------------------------------------


gwas = fread(sumstats.file)

# first remove any snps not passing threshold
# and only need `#CHR`, `BP`, `P`
gwas = gwas[P < pval.th][, .(`#CHR`, BP, P)] %>% unique
gwas[, prank := rank(P, ties.method = "first")]
gwas = gwas[order(prank)]

print(paste0("Found ", nrow(gwas), " SNPs passing threshold."))

print(glue("Start finding GWAS hit loci with pval < {pval.th}."))
print(glue("{nrow(gwas)} SNPs passing threshold to be processed."))

# make GRanges
gwas.gr = makeGRangesFromDataFrame(
    gwas,
    keep.extra.columns = T,
    ignore.strand = T,
    seqnames.field = "#CHR",
    start.field = "BP",
    end.field = "BP")


# Find Loci --------------------------------------------------------------------


# initialize
max.loop = length(gwas.gr) # prevent infinite loop
locusSize = as.integer(1e6) 
loci = list()
counter = 1
# Loop over all the considered GWAS SNPs ton construct loci
while (length(gwas.gr) > 0 & counter < max.loop) {
# while (length(gwas.gr) > 0 & counter < 10) {
    res = findGwasHitLocus(gwas.gr, locusWidth = locusSize)
    if (!is.null(res$HitLocus)) {
        loci = c(loci, res$HitLocus)
    }
    gwas.gr = res$RemainedSnps
    counter = counter + 1
    if (counter %% 100 == 0) {
        print(glue("Processed {counter} Loci."))
    }
}
print(glue("Processed {counter} Loci."))

# convert list of granges to a single grange object and sort
loci = map_dfr(loci, as.data.frame) %>% 
    makeGRangesFromDataFrame(keep.extra.columns = T) %>% 
    sort

# reduce overlapping loci
print(glue("Merge overlapping loci."))
loci = GenomicRanges::reduce(loci)
print(glue("Found {length(loci)} loci. Summary of Loci width:\n"))
print(summary(width(loci)))


# export out as bed file
print(glue("Exporting loci to {outHits.file}."))
rt = rtracklayer::export.bed(loci, outHits.file)


