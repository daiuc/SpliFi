# prepare input data for dge (GTEs), namely
# 1. raw counts
# 2. column data

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
    stop("Usage: Rscript prepare_GTEx_dge.R <tissue1_count_file> <tissue2_count_file> <number_of_samples> <outprefix>")
}

tissue1_count_file <- args[1]
tissue2_count_file <- args[2]
number_of_samples <- as.numeric(args[3])
outprefix <- args[4]

suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(glue))

print(glue("input arguments: 
    tissue1_count_file: {tissue1_count_file}
    tissue2_count_file: {tissue2_count_file}
    number_of_samples: {number_of_samples}
    outprefix: {outprefix}"))

t1_counts <- fread(tissue1_count_file, header = TRUE, sep = "\t")
t2_counts <- fread(tissue2_count_file, header = TRUE, sep = "\t")
tissue1 <- str_split(tissue1_count_file, "/") %>%
    unlist %>%
    tail(1) %>%
    str_remove("_gene_reads.tsv.gz")

tissue2 <- str_split(tissue2_count_file, "/") %>%
    unlist %>%
    tail(1) %>%
    str_remove("_gene_reads.tsv.gz")


print(glue("Randomly select {number_of_samples} samples from {tissue1} and {tissue2}"))

# stop if number_of_samples is greater than ncol - 2 for either t1_counts or t2_counts
if (number_of_samples > (ncol(t1_counts) - 2) | number_of_samples > (ncol(t2_counts) - 2)) {
    stop(glue("Exiting, trying to subset more samples than available in either {tissue1} or {tissue2}!\n"))
}

t1_data_cols <- colnames(t1_counts)[3:ncol(t1_counts)]
t1_selected_cols <- sample(t1_data_cols, number_of_samples)
t2_data_cols <- colnames(t2_counts)[3:ncol(t2_counts)]
t2_selected_cols <- sample(t2_data_cols, number_of_samples)

print(glue("Selected columns for {tissue1}: {paste(t1_selected_cols, collapse = ', ')}"))
print("\n\n")
print(glue("Selected columns for {tissue2}: {paste(t2_selected_cols, collapse = ', ')}"))

t1_selected_counts <- t1_counts[, c('Name', 'Description', t1_selected_cols), with = FALSE]
t2_selected_counts <- t2_counts[, c('Name', 'Description', t2_selected_cols), with = FALSE]

# combine the selected counts
combined_counts <- cbind(t1_selected_counts, t2_selected_counts[, -c('Name', 'Description')])

# col data
coldata <- data.table(
    sample_id = c(t1_selected_cols, t2_selected_cols),
    tissue = c(rep(tissue1, number_of_samples), rep(tissue2, number_of_samples))
)

out_count_file <- glue("{outprefix}/{tissue2}_v_{tissue1}_counts.tsv")
out_coldata_file <- glue("{outprefix}/{tissue2}_v_{tissue1}_coldata.tsv")

fwrite(combined_counts, out_count_file, sep = "\t")
fwrite(coldata, out_coldata_file, sep = "\t")

print(glue("Done! Output files: 
counts written to: {out_count_file}
Column data written to: {out_coldata_file}"))





