# This is copied directly from Ben's PrepareGTExPhenotypes.R script

library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(edgeR)
library(RNOmni)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
    stop("Usage: Rscript PrepareGTExPhenotypes.R <data.path> <genes.path> <CPM.Bed_Out> <QQnorm.Bed_Out>")
}

data.path <- args[1]
genes.path <- args[2]
CPM.Bed_Out <- args[3]
QQnorm.Bed_Out <- args[4]



# data.path <- 'GTEx/data/gene_reads_2017-06-05_v8_adipose_visceral_omentum.gct.gz'

# dat <- read_tsv(data.path, skip = 2)
dat <- fread(data.path)

genes <- read_tsv(genes.path,
    col_names = c("Chr", "Start", "End", "Geneid", "score", "strand")
) %>%
    mutate(gene = str_extract(Geneid, "^[^.]+"))

dat <- dat %>%
    mutate(gene = str_extract(Name, "^[^.]+")) %>%
    filter(gene %in% genes$gene) %>%
    select(-Description, -Name) %>%
    column_to_rownames(., "gene") %>%
    filter(rowSums(. != 0) > 0)

colnames(dat) <- str_extract(colnames(dat), "Name|Description|GTEX\\-[\\w\\d]+")

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
