#!/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(glue))
suppressMessages(library(data.table))


splitFrac <- function(x) {
    x  <- str_split(x, "/")
    return(as.integer(map(x, ~.x[1])))
}

processDT <- function(dt) {
    chroms <- dt$chrom
    datacols <- names(dt)[-c(1)]
    dt <- dt[, -c('chrom')][, map(.SD, ~splitFrac(.x))]
    dt$chrom <- chroms
    dt$intron_type = str_extract(dt$chrom, '[FN]')
    dt$clu = str_extract(dt$chrom, 'clu_\\d+')
    dt[, clu_type := paste(unique(intron_type), collapse = ","), by = clu]
    cols <- c('chrom', 'clu', 'clu_type', 'intron_type', datacols)

    return(dt[, ..cols])

}


args <- commandArgs(trailingOnly = TRUE)

tissue  <- args[1]

count_file  <- glue('/project2/yangili1/cdai/SpliFi/code/results/pheno/noisy/GTEx/{tissue}/wConst_perind.constcounts.noise_by_intron.gz')
out_file  <- glue('data-for-2023-12-05-Notebook-{tissue}-allIntrons.rds')
counts <- fread(count_file)
counts <- processDT(counts)

saveRDS(counts, out_file)
