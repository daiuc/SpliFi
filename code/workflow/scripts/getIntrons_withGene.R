# get intron targets to plot for sashimi plots


suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(glue))

if (!interactive()) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) < 2) {
        stop("Run as: Rscript getIntrons.R <data.rds> <contrast> <outpath>> <minDeltaPSI:default 0.3> <minPvalue:default 1e-3> <minL2FC:default 0.58> <genelist.txt>")
    }
} else {
    args <- c(
        "/project/yangili1/cdai/SpliFi/data/ds_v_dge/Brain-Cortex_v_Muscle-Skeletal_data.rds",
        "Brain-Cortex_v_Muscle-Skeletal",
        "./out-test.txt",
        "0.2", # min deltaPSI for ds
        "1e-3", # fdr for ds and dge
        "0.58", # min log2FoldChange for dge
        "test.genes.txt" # file with gene names to select from one gene name per line, no header
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
select_genes <- read_lines(args[7]) # gene names to select/keep

print(glue("Get introns for sashimi plots for {contrast} with:\nFDR cutoff {FDR},\nmin delta PSI cutoff {minDeltaPsi},\nmin dge l2fc > {minL2FC}.\nSelect introns from {length(select_genes)} genes defined in {args[7]}"))

ds_keep_cols <- c("intron", "cluster", "itype", "ctype", "gene_id", "gene_name", "deltapsi", "p.adjust")
plotDt <- inner_join(
    ds[ctype == "PR,UP" & itype == "UP" & gene_name %in% select_genes, ..ds_keep_cols],
    dge[, .(gene_id, log2FoldChange, `p.adjust` = padj)],
    by = "gene_id", suffix = c("_ds", "_dge")
)

plotDt <- plotDt[
    `p.adjust_ds` < FDR &
        `p.adjust_dge` < FDR &
        abs(deltapsi) > minDeltaPsi &
        abs(log2FoldChange) > minL2FC &
        deltapsi * log2FoldChange < 0,
] %>%
    .[order(-abs(deltapsi)), .(intron,
        deltapsi = round(deltapsi, 5),
        l2fc = round(log2FoldChange, 5), gene_id, itype,
        ctype, `p.adjust_ds`, `p.adjust_dge`
    )]

# write intron to lines
if (nrow(plotDt) == 0) {
    stop(glue("No introns found for defined selection criteria."))
} else {
    print(glue("Writing {nrow(plotDt)} introns to {outFile}"))
    plotDt[, intron] %>%
        unique() %>%
        writeLines(outFile)
}
# fwrite(plotDt, outFile, sep = '\t', quote = FALSE, row.names = FALSE)
