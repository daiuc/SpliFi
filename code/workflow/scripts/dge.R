# differential gene expression analysis using DESeq2
# currently specifically designed for GTEx

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
    stop("Usage: Rscript dge.R <raw_counts_file> <column_data_file> <outprefix> [min_reads (default: 10)] [min_samples: (default: 5)]")
}
raw_counts_file <- args[1]
column_data_file <- args[2]
outprefix <- args[3]
min_reads <- ifelse(length(args) > 3, as.numeric(args[4]), 10)
min_samples <- ifelse(length(args) > 4, as.numeric(args[5]), 5)


suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(DESeq2))
suppressMessages(library(glue))

out_df_file <- paste0(outprefix, "_dge_genes.tsv")
out_description_file <- paste0(outprefix, "_dge_description.txt")


print(glue("input arguments: 
    raw_counts_file: {raw_counts_file}
    column_data_file: {column_data_file}
    min_reads: {min_reads}
    min_samples: {min_samples}"))

print(glue("output files: 
    Result table: {out_df_file}
    Result description: {out_description_file}"))


# load input data
# NOTE: two input data files are required, both in .tsv format:
# 1. gene expression raw counts, first column is gene id, subseqent columns are read counts for each sample
# gene expression data has columns: Name, Description, and then sample ids
# 2. column data, first column is sample id, all rows of the first column match exactly 
# the 2nd... nth column of the gene expression raw counts file
# 3. second column in column data is tissue type.
# 4. both raw counts and column data files' first column will be turned into row-names

raw_counts <- fread(raw_counts_file, header = TRUE, sep = "\t")
raw_counts <- raw_counts[, -c('Description')]
raw_counts <- column_to_rownames(raw_counts, var = "Name")

column_data <- fread(column_data_file, header = TRUE, sep = "\t")
column_data <- column_to_rownames(column_data, var = "sample_id")
tissues <- column_data[, 'tissue'] %>% unique
print("Tissues: paste(tissues, collapse = ', ')")


# construct a deseq dataset object
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = column_data,
                              design = ~ tissue)

print(glue("dds object constructed with {nrow(dds)} genes and {ncol(dds)} samples"))
print(glue("contrast: {tissues[2]} (numerator) vs. {tissues[1]} (denominator)"))
dds$tissue <- factor(dds$tissue, levels = tissues)

# pre-filtering
print(glue("pre-filtering: removing genes with less than {min_reads} reads in at least {min_samples} samples"))
keep  <- rowSums(counts(dds) >= min_reads) >= min_samples
dds <- dds[keep,]
print(glue("after pre-filtering, {nrow(dds)} genes remain"))

# run DESeq2
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(-abs(res$log2FoldChange), res$pvalue), ]

# write results to file (tsv)
print(glue("writing {nrow(res)} genes to {out_df_file}. Column descriptions in {out_description_file}"))
res.dt <- as.data.table(res, keep.rownames = 'gene_id')
fwrite(res.dt, out_df_file, sep = "\t")


# write description to file
write_lines(res@elementMetadata@listData, out_description_file)


