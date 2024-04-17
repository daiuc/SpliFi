# get intron targets to plot for sashimi plots


suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(glue))

if (!interactive()) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) < 2) {
        stop('Run as: Rscript getIntrons.R <data.rds> <contrast> <outpath>> <minDeltaPSI:default 0.3> <minPvalue:default 1e-3> <minL2FC:default 0.58>.')
    }
} else {
  args <- c(
            '/project2/yangili1/cdai/SpliFi/data/ds_v_dge/Brain-Cortex_v_Muscle-Skeletal_data.rds',
            'Brain-Cortex_v_Muscle-Skeletal',
            './out-test.txt',
            '0.2', # min deltaPSI for ds
            '1e-3', # fdr for ds and dge
            '0.58' # min log2FoldChange for dge
            )
}

# input data
data <- readRDS(args[1])
ds <- data$ds
dge <- data$dge
contrast <- args[2]
outFile <- args[3]
minDeltaPsi <- if (is.na(args[4])) 0.2 else as.numeric(args[4])
FDR <- if (is.na(args[5])) 1e-3 else as.numeric(args[5])
minL2FC <- if (is.na(args[6])) 0.58 else as.numeric(args[6])

print(glue('Get introns for sashimi plots for {contrast} with FDR cutoff {FDR},  min delta PSI cutoff {minDeltaPsi}, and min dge l2fc > {minL2FC}'))

ds_keep_cols <- c('intron', 'cluster', 'itype', 'ctype', 'gene_id', 'deltapsi', 'p.adjust')
plotDt <- inner_join(
    x = ds[ctype == 'PR,UP' & itype == 'UP', ..ds_keep_cols],
    y = dge[, .(gene_id, log2FoldChange, `p.adjust` = padj)],
    by = 'gene_id', suffix = c('_ds', '_dge')
    ) 

plotDt <- plotDt[
        `p.adjust_ds` < FDR &
        `p.adjust_dge` < FDR &
        abs(deltapsi) > minDeltaPsi &
        abs(log2FoldChange) > minL2FC &
        deltapsi * log2FoldChange < 0,
    ] %>% 
    .[order(-abs(deltapsi)), .(intron, deltapsi = round(deltapsi, 5), 
        l2fc = round(log2FoldChange, 5), gene_id, itype, 
        ctype, `p.adjust_ds`, `p.adjust_dge`)]

# write intron to lines
print(glue('Writing {nrow(plotDt)} introns to {outFile}'))
plotDt[, intron] %>% unique() %>% writeLines(outFile)
# fwrite(plotDt, outFile, sep = '\t', quote = FALSE, row.names = FALSE)





