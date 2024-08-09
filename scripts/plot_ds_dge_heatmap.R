# this script is intended to be run interactively only.


if (! interactive()) {
  stop("This script is intended to be run interactively only.")
}

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(glue))

suppressMessages(library(cowplot))
theme_set(theme_cowplot(font_size = 14))

library(httpgd)

setwd("~/cdai/SpliFi")
getwd()

hgd(host = '10.50.250.58', port = 9999)


#--- functions


GetHeatmapMatrix <- function(data, introns, contrast) {
  data$intron2 <- str_replace_all(data$intron, ":clu_\\d+_", "")
  tissues <- str_split(contrast, "_v_") %>% unlist
  keep1 <- c('intron2', tissues[1]) 
  keep2 <- c('intron2', tissues[2]) 
  out <- list(
    data[intron2 %in% introns, keep1, with = FALSE],
    data[intron2 %in% introns, keep2, with = FALSE]
  )
  names(out) <- tissues
  
  return(out)

}

GetTopUpIntrons <- function(data, FDR_ds, FDR_dge, dPSI) {
  dge <- data$dge
  ds <- data$ds

  # if a cluster has multiple UP introns, only select the best 1
  ds <- ds[itype == 'UP' & ctype == 'PR,UP'][, rk := rank(-abs(deltapsi), ties.method = "first"), by = cluster][rk ==1][, rk := NULL][]
  ds <- ds[`p.adjust` < FDR_ds & abs(deltapsi) > dPSI,]

  dge <- dge[padj < FDR_dge,]
  
  ds_excl_cols <- c('itype', 'ctype', 'df', 'p', 'p.adjust','logef', 'loglr', 'status')
  dge_excl_cols <- c('baseMean', 'lfcSE', 'stat', 'pvalue', 'padj')
  chosen <- inner_join(
      x = ds[, -ds_excl_cols, with = FALSE],
      y = dge[, -dge_excl_cols, with = FALSE],
      by = "gene_id",
      suffix = c("_ds", "_dge")
    ) %>%
    .[deltapsi * log2FoldChange < 0, ]

  return(chosen)
}


#--- read in prepared data
contrast_ls <- dir('~/cdai/SpliFi/data/ds_v_dge', '*.rds') %>% str_remove_all('_data\\.rds')
tissues <- str_split(contrast_ls, "_v_") %>% unlist %>% unique

data_f <- glue("~/cdai/SpliFi/data/ds_v_dge/{contrast_ls}_data.rds")
names(data_f) <- contrast_ls
data <- map(data_f, readRDS)

#-- play

length(contrast_ls)
head(tissues)




