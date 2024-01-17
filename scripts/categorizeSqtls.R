#!/bin/env Rscript
suppressMessages(library(tidyverse))
suppressMessages(library(glue))
suppressMessages(library(data.table))


inFileList <- commandArgs(trailingOnly = TRUE)[1]
fdr <- commandArgs(trailingOnly = TRUE)[2]
fdr <- as.numeric(fdr)
outFile <- commandArgs(trailingOnly = TRUE)[3]


# categorize sQTL to u-sQTL or p-sQTL
cat_sqtl <- function(dt) {
    # each dt is a qtl_sig subsetted by a cluster id based on columns: clu_cl, pval_nom
    if (dt[1, clu_cl] == "N") {
        # 1. if clu_cl == "N", then u-sQTL
        dt$qtl_cl <- "u-sQTL"
        dt <- dt[order(-abs(slope))][1] # select one with biggest nominal effect
    } else if (dt[1, clu_cl] == "F") {
        # 2. if clu_cl == "F", then p-sQTL
        dt$qtl_cl <- "p-sQTL"
        dt <- dt[order(-abs(slope))][1] # select one with biggest nominal effect
    } else {
        # 3. if clu_cl == "N,F", then label as u-sQTL
        dt$qtl_cl <- "u-sQTL"
        dt <- dt[junc_cl == "N"][order(-abs(slope))][1] # select biggest unprod nominal effect
    }

    return(dt)
}





# Read file path from inFileList
cat("Reading files from:\n", inFileList, "\n")
qtl_files <- readLines(inFileList)
names(qtl_files) <- str_extract(qtl_files, "chr\\d{1,2}")
# sort files by name use natural sort
qtl_files <- naturalsort::naturalsort(qtl_files)

# read in files
cat("Reading files:\n", paste(qtl_files, collapse = "\n"), "\n")
qtl_dt <- map_dfr(qtl_files, fread)

# add cluster_id
cat("Adding cluster id\n")
qtl_dt[, `:=`(
              # cluster id
              clu_id = str_extract(phenotype_id, "clu_\\d+"),
              # junction type (either N for unproductive or F for productive)
              junc_cl = str_sub(phenotype_id, -1))]


# classify cluster
cat("Classifying clusters\n")
qtl_dt[, clu_cl := paste(sort(unique(junc_cl)), collapse = ","), by = clu_id]

cat(glue("Using min FDR: {fdr}\n\n"))
# store significant sQTLs
# 1. remove clusters with only unproductive junctions
qtl_sig <- qtl_dt[ clu_cl != "N"]
# cluster ids with at least 1 significant splice junction
sqtl_clu_ids <- qtl_sig[, .(N_sig = sum(q < fdr)), by = clu_id][N_sig > 0, clu_id]
# keep only significant clusters
qtl_sig <- qtl_sig[clu_id %in% sqtl_clu_ids]
# keep only significant introns
qtl_sig <- qtl_sig[q < fdr]

# re-lable clu_cl based on significant junctions only
cat("Re-labelling clusters\n")
qtl_sig[, clu_cl := paste(sort(unique(junc_cl)), collapse = ","), by = clu_id]

# categorize significant junctions to u-sQTL or p-sQTL, 1 per cluster
cat("Categorizing sQTLs\n")
qtl_sig <- split(qtl_sig, by = "clu_id")  %>% map_dfr(cat_sqtl)


# write to file
cat("Writing to file:\n", outFile, "\n")
fwrite(qtl_sig, outFile, sep = "\t")

