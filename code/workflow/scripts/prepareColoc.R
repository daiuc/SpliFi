# Prepare coloc input files
# NOTE: eqtl phenotype is based on count table which uses gene_ids from gencode.v26.GRCh38.genes.gtf
#       Make sure when looking up intron's gene_id, use the same gtf file.

if (interactive()) {
    args <- scan(
        text = "results/eqtl/GTEx/Liver/perm/chr1.txt.addQval.txt.gz
                         results/qtl/noisy/GTEx/Liver/separateNoise/cis_100000/perm/chr1.addQval.txt.gz
                         0.1
                         /project/yangili1/cdai/annotations/hg38/gencode.v26.GRCh38.genes.csv
                         results/eqtl/GTEx/Liver/qqnorm.sorted.chr1.bed.gz
                         results/pheno/noisy/GTEx/Liver/separateNoise/leafcutter.qqnorm_chr1.gz
                         results/coloc/sqtl-eqtl/GTEx/Liver/chr1.eqtl_pheno.bed.gz
                         results/coloc/sqtl-eqtl/GTEx/Liver/chr1.sqtl_pheno.bed.gz
                         results/coloc/sqtl-eqtl/GTEx/Liver/chr1.eqtl_sqtl_ids.txt
                        ",
        what = character()
    )
    eqtl_perm_f <- args[1]
    sqtl_perm_f <- args[2]
    FDR <- as.numeric(args[3])
    genes_f <- args[4]
    eqtl_pheno_f <- args[5]
    sqtl_pheno_f <- args[6]
    eqtl_pheno_o <- args[7]
    sqtl_pheno_o <- args[8]
    coloc_pheno_ids_f <- args[9]
} else {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) < 2) {
        stop("Two arguments are required")
    }
    eqtl_perm_f <- args[1]
    sqtl_perm_f <- args[2]
    FDR <- as.numeric(args[3])
    genes_f <- args[4]
    eqtl_pheno_f <- args[5]
    sqtl_pheno_f <- args[6]
    eqtl_pheno_o <- args[7]
    sqtl_pheno_o <- args[8]
    coloc_pheno_ids_f <- args[9]
}

suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(GenomicRanges))
library(glue)

#--- load perm pass qtls

print(glue("Load perm pass qtls: \neqtl: {eqtl_perm_f}, \nsqtl: {sqtl_perm_f}"))
eqtl_perm <- fread(eqtl_perm_f)
sqtl_perm <- fread(sqtl_perm_f)


#--- filter qtls
print(glue("Filter qtls with FDR < {FDR}"))
eqtl_perm <- eqtl_perm[q < FDR]
sqtl_perm <- sqtl_perm[q < FDR & str_detect(phenotype_id, "PR|UP|NE")]

print(glue("Number of QTLs after filtering: \neQTL: {nrow(eqtl_perm)}, \nsQTL: {nrow(sqtl_perm)}"))

#--- load gene annotation
print(glue("Load gene annotation: {genes_f}"))
genes <- fread(genes_f)
genes <- genes[feature == "gene", .(seqname, start, end, gene_id, score, strand)] %>% unique()
print(glue("{nrow(genes)} genes are loaded."))

#--- construct intron_id, gene_id lookup table based on distance
# gene must complete enclose intron
sqtl.gr <- sqtl_perm[, .(seqname = phenotype_chr, start = phenotype_start, end = phenotype_end, strand = phenotype_strand, intron_id = phenotype_id)] %>%
    makeGRangesFromDataFrame(keep.extra.columns = T, starts.in.df.are.0based = F)
genes.gr <- makeGRangesFromDataFrame(genes, keep.extra.columns = T, starts.in.df.are.0based = F)

intron_to_gene <- distanceToNearest(sqtl.gr, genes.gr, select = "all", ignore.strand = F)
rowindex <- which(intron_to_gene@elementMetadata$distance == 0)

df1 <- sqtl.gr[queryHits(intron_to_gene[rowindex])] %>%
    as.data.frame() %>%
    as.data.table() %>%
    .[, .(c1 = seqnames, s1 = start, e1 = end, st1 = strand, intron_id = intron_id)]
df2 <- genes.gr[subjectHits(intron_to_gene[rowindex])] %>%
    as.data.frame() %>%
    as.data.table() %>%
    .[, .(c2 = seqnames, s2 = start, e2 = end, st2 = strand, gene_id = gene_id)]

lookup <- cbind(df1, df2) %>%
    .[s1 >= s2 & e1 <= e2 & st1 == st2, .(intron_id, gene_id)]

#--- use lookup table to get gene_id for each intron
sqtl_perm <- inner_join(sqtl_perm, lookup, by = c("phenotype_id" = "intron_id"))

#--- select eqtl and sqtl phenotypes

select_phenos <- sqtl_perm[, .(phenotype_id, gene_id)] %>% unique()
print(glue("Selected {nrow(select_phenos)} sQTL-eQTL pairs for coloc,"))
print(glue("with {uniqueN(select_phenos$gene_id)} unique genes and {uniqueN(select_phenos$phenotype_id)} unique introns."))


# remove gene_id suffix of .1, et.c
print(glue("Save selected phenotype eQTL and sQTL phenotype ids to {coloc_pheno_ids_f}"))
select_phenos[, gene_id := str_replace_all(gene_id, "\\..+$", "")]
fwrite(select_phenos, coloc_pheno_ids_f, sep = "\t")

#--- write selected splicing and expression phenotypes to table for nominal pass qtltools
# use selected phenotype_id to filter qqnormalized leafcutter output
sqtl_pheno <- fread(sqtl_pheno_f)
sqtl_pheno <- sqtl_pheno[pid %in% select_phenos$phenotype_id]

eqtl_pheno <- fread(eqtl_pheno_f)
eqtl_pheno$pid <- str_replace_all(eqtl_pheno$pid, "\\..+$", "")
eqtl_pheno <- eqtl_pheno[pid %in% select_phenos$gene_id]
# eqtl_pheno$pid <- eqtl_pheno$gid


# make sure properly sorted before write out
chroms <- c(glue("chr{1:22}"), "chrX", "chrY")
sqtl_pheno[, `#Chr` := factor(`#Chr`, levels = chroms)]
sqtl_pheno <- sqtl_pheno[order(`#Chr`, start, end)]
eqtl_pheno[, `#Chr` := factor(`#Chr`, levels = chroms)]
eqtl_pheno <- eqtl_pheno[order(`#Chr`, start, end)]

# use these two output pheno bed to run qtltools for nominal pass
print(glue("Save selected eQTL and sQTL phenotypes to {eqtl_pheno_o} and {sqtl_pheno_o} for qtltools nominal pass."))
fwrite(sqtl_pheno, sqtl_pheno_o, sep = "\t")
fwrite(eqtl_pheno, eqtl_pheno_o, sep = "\t")

print("Done.")


### scratch

# dim(eqtl_perm)
# dim(sqtl_perm)
#
# eqtl_perm[1:2, ]
# sqtl_perm[1:2, ]
#
# head(genes)
# dim(genes)
#
# sqtl.gr
# genes.gr
#
# sqtl_pheno %>% dim()
# sqtl_pheno[1:2, 1:7]
#
#
# eqtl_pheno %>% dim()
# eqtl_pheno[1:2, 1:7]
#
# select_phenos$gene_id %>% uniqueN()
# select_phenos
