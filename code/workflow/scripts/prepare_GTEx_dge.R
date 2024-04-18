# prepare input data for dge (GTEs), namely
# 1. raw counts
# 2. column data

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
    stop("Usage: Rscript prepare_GTEx_dge.R <tissue1_count_file> <tissue2_count_file> <ds_sample_file> <outprefix>")
}

tissue1_count_file <- args[1]
tissue2_count_file <- args[2]
ds_sample_file <- args[3]
outprefix <- args[4]

if (interactive()) {
    tissue2_count_file <- 'resources/GTEx/expression/Brain-Cortex_gene_reads.tsv.gz'
    tissue1_count_file <- 'resources/GTEx/expression/Muscle-Skeletal_gene_reads.tsv.gz'
    ds_sample_file <- 'results/ds/GTEx/Brain-Cortex_v_Muscle-Skeletal/ds_sample_group.txt'
    outprefix <- 'test'
}

suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(glue))

print(glue("input arguments: 
    tissue1_count_file: {tissue1_count_file}
    tissue2_count_file: {tissue2_count_file}
    sample_file_from_DS: {ds_sample_file}
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

ds_samples <- fread(ds_sample_file, header = FALSE, sep = " ", col.names = c('sample_id', 'group'))
ds_samples$ind_id <-  str_extract(ds_samples$sample_id, "GTEX-[\\d\\w]+")

# check group/tissue names are the same
if (length(intersect(c(tissue1, tissue2), ds_samples$group)) != 2) {
    stop(glue("Exiting, tissue names in {ds_sample_file} do not match the input files!\n"))
}

t1_data_cols <- data.table(cols = colnames(t1_counts)[3:ncol(t1_counts)])
t1_data_cols[, ind_id := str_extract(cols, "GTEX-[\\d\\w]+")]
t1_selected_cols <- inner_join(t1_data_cols, ds_samples[group == tissue1], by = 'ind_id')$cols

t2_data_cols <- data.table(cols = colnames(t2_counts)[3:ncol(t2_counts)])
t2_data_cols[, ind_id := str_extract(cols, "GTEX-[\\d\\w]+")]
t2_selected_cols <- inner_join(t2_data_cols, ds_samples[group == tissue2], by = 'ind_id')$cols

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
    tissue = c(rep(tissue1, length(t1_selected_cols)), rep(tissue2, length(t2_selected_cols))))

out_count_file <- glue("{outprefix}/{tissue2}_v_{tissue1}_counts.tsv")
out_coldata_file <- glue("{outprefix}/{tissue2}_v_{tissue1}_coldata.tsv")

fwrite(combined_counts, out_count_file, sep = "\t")
fwrite(coldata, out_coldata_file, sep = "\t")

print(glue("Done! Output files: 
counts written to: {out_count_file}
Column data written to: {out_coldata_file}"))





