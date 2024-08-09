# This is copied directly from Ben's PrepareGTExPhenotypes.R script

library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(edgeR)
library(RNOmni)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
    stop("Usage: Rscript PrepareGTExPhenotypes.R <data.path> <genes.path> <samples.path> <CPM.Bed_Out> <QQnorm.Bed_Out>")
}

data.path <- args[1]
genes.path <- args[2]
samples.path <- args[3]
CPM.Bed_Out <- args[4]
QQnorm.Bed_Out <- args[5]

if (interactive()) {
    args <- scan(text = "Rscript workflow/scripts/fromBen_PrepareGTExEQTLPheno.R resources/GTEx/expression/Brain-Amygdala_gene_reads.tsv.gz resources/GTEx/ExpressedGeneList.txt results/pheno/noisy/GTEx/Brain-Amygdala/separateNoise/leafcutter_names.txt results/eqtl/GTEx/Brain-Amygdala/cpm.bed.gz results/eqtl/GTEx/Brain-Amygdala/qqnorm.bed.gz", what = character())
    data.path <- args[3]
    genes.path <- args[4]
    samples.path <- args[5]
}


dat <- fread(data.path)

genes <- read_tsv(genes.path,
    col_names = c("Chr", "Start", "End", "Geneid", "score", "strand")
) %>%
    mutate(gene = str_extract(Geneid, "^[^.]+"))

samples <- read_lines(samples.path)

dat <- dat %>%
    mutate(gene = str_extract(Name, "^[^.]+")) %>%
    filter(gene %in% genes$gene) %>%
    select(-Description, -Name) %>%
    column_to_rownames(., "gene") %>%
    filter(rowSums(. != 0) > 0)

colnames(dat) <- str_extract(colnames(dat), "Name|Description|GTEX\\-[\\w\\d]+")
kept_samples <- intersect(samples, colnames(dat))

dat <- dat %>%
    select(all_of(kept_samples))
keepRows <- rowSums(dat) >= 10
dat <- dat[keepRows, ]

dat.cpm <- dat %>%
    cpm(log = T, prior.count = 0.1)

genes <- genes %>% filter(gene %in% rownames(dat.cpm))

dat.standardized <- dat.cpm %>%
    t() %>%
    scale() %>%
    t()
dat.qqnormed <- apply(dat.standardized, 2, RankNorm)

out.cpm <- genes %>%
    inner_join((dat.cpm %>% as.data.frame() %>% rownames_to_column("gene"))) %>%
    mutate(start = as.numeric(Start)) %>%
    mutate(across(where(is.numeric), round, 5)) %>%
    dplyr::select(`#Chr` = Chr, start, end = End, pid = Geneid, gid = Geneid, strand = strand, everything(), -Start, -score, -gene) %>%
    arrange(`#Chr`, start)

out.qqnormed <- genes %>%
    inner_join((dat.qqnormed %>% as.data.frame() %>% rownames_to_column("gene"))) %>%
    mutate(start = as.numeric(Start)) %>%
    mutate(across(where(is.numeric), round, 5)) %>%
    dplyr::select(`#Chr` = Chr, start, end = End, pid = Geneid, gid = Geneid, strand = strand, everything(), -Start, -score, -gene) %>%
    arrange(`#Chr`, start)

write_tsv(out.cpm, CPM.Bed_Out)
write_tsv(out.qqnormed, QQnorm.Bed_Out)
