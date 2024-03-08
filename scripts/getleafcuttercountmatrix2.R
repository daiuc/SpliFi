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
    dt$intron_type = str_sub(dt$chrom, -2, -1)
    dt$clu = str_extract(dt$chrom, 'clu_\\d+')
    dt[, clu_type := paste(sort(unique(intron_type)), collapse = ","), by = clu]
    cols <- c('chrom', 'clu', 'clu_type', 'intron_type', datacols)

    return(dt[, ..cols])

}


args <- commandArgs(trailingOnly = TRUE)

inFile  <- args[[1]]
outFile  <- args[[2]]

print(inFile)
print(outFile)

counts <- fread(inFile)
counts <- processDT(counts)

saveRDS(counts, outFile)
print(glue("done at:{Sys.time()}"))

