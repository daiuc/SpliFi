# This script is used to prepare data needed for analyses and plots
# downostream of teh DGE analysis. 

# NOTE: Mainly it produces these data: 
# For differential splicing: 
# 1. A table of effect sizes, with intron cluster id column added
# 2. A table of intron id, cluster id, intron type, and cluster type
# 3. A table of cluser level significance. 
# For differential gene expression:


if (interactive()) {
  contrast <- "Liver_v_Lung"
  gtf.f <- "/project2/yangili1/cdai/annotations/hg38/gencode.v26.GRCh38.genes.csv"
  dsPrefix <- "results/ds/GTEx"
  dgePrefix <- "results/dge/GTEx"

} else {
  args = commandArgs(trailingOnly = TRUE)
  contrast <- args[1]
  gtf.f <- args[2]
  dsPrefix <- args[3]
  dgePrefix <- args[4]
  outPrefix <- args[5]

  if (length(args) < 5) {
    stop("Usage: Rscript prepDGE_DS_AnalysesData.R <contrast: eg. Liver_v_Lung> <gencode v26 gtf csv file> <ds prefix> <dge prefix> <output prefix>")
  }

  print(glue::glue("contrasts: {contrast}
    gtf: {gtf.f}
    dsPrefix: {dsPrefix}
    dgePrefix: {dgePrefix}
    outPrefix: {outPrefix}"))
}


suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(data.table))
suppressMessages(library(bedtoolsr))
suppressMessages(library(GenomicRanges))

# ----------------- Helper functions -----------------
# get intron labels
labelIntron <- function(df) {
    df <- df[, 
          # get intron id (removing label)
        .(intron = str_sub(chrom, 1, -4),
          # extract cluster id
          cluster = str_split(chrom, ":") %>% map_chr(~paste(.x[1], .x[4], sep=":")),
          # extract intron type (label)
          itype = str_sub(chrom, -2, -1)
    )]
    # get cluster type
    df <- df[, ctype := paste(sort(unique(itype)), sep="", collapse=","), by = cluster][]

    return(df)
}

# ----------------- Load Gencode ----------------- load gtf for gene_id, gene_name gtf <- fread(gtf)
gtf <- fread(gtf.f) %>% 
  .[feature == 'gene' & gene_type == 'protein_coding',
    .(seqname, start, end, gene_name, gene_id, strand)] %>%
  unique()

gtf <- makeGRangesFromDataFrame(gtf,
    keep.extra.columns = TRUE,
    ignore.strand = FALSE
)

# ----------------- Differential Splicing data -----------------

ds.pval <- fread(glue::glue("{dsPrefix}/{contrast}/ds_cluster_significance.txt"))
ds.beta <- fread(glue::glue("{dsPrefix}/{contrast}/ds_effect_sizes.txt"))
# only grab the first column, which is intron id with labels
ds.intron <- fread(glue::glue("{dsPrefix}/{contrast}/ds_perind_numers.counts.noise_by_intron.gz"))[, 1] 

# add cluster id and intron type
ds.intron <- labelIntron(ds.intron)

# use intron coordinates to get gene_name and gene_id
coords <- str_split(ds.intron$intron, ":", simplify = TRUE)  %>% as.data.table()
setnames(coords, c("seqname", "start", "end", "cluster"))
coords[, strand := str_sub(cluster, -1, -1)]
coords[, cluster := NULL]
ds.intron <- cbind(ds.intron, coords)
ds.intron <- makeGRangesFromDataFrame(
  ds.intron, 
  keep.extra.columns = TRUE,
  ignore.strand = FALSE,
  starts.in.df.are.0based = TRUE)

# use overlap to get gene_name and gene_id that intron belongs to
olaps <- findOverlaps(ds.intron, gtf, type="within", select="all", ignore.strand=FALSE)

an.intron <- ds.intron[olaps@from] # annotated introns
mcols(an.intron) <- cbind(mcols(gtf[olaps@to]), mcols(an.intron))

# remove introns that are mapped to multiple genes
an.intron <- as.data.table(an.intron)[
  , .(gene_name, gene_id, rk = frank(gene_name)), 
  by = .(seqnames, start, end, strand, intron, cluster, itype, ctype)
  ][, .(gene_name, gene_id, maxrk = max(rk)), 
    by = .(seqnames, start, end, strand, intron, cluster, itype, ctype)
  ][maxrk == 1, -c("maxrk")]

# merge with ds.beta to bring in ds effect size (intron level)
an.intron <- inner_join(an.intron[, -c("seqnames", "start", "end", "strand")], ds.beta, by = c("intron"))

# merge with ds.pval to bring in ds p-value (cluster level)
an.intron <- inner_join(an.intron, ds.pval, by = c("cluster"))


# ----------------- Differential Expression data -----------------

# read in vectors from  a text file
dge.dsc <- read_lines(glue::glue("{dgePrefix}/{contrast}_dge_description.txt"))
dge.dsc <- eval(parse(text = dge.dsc[2]))

# dge results
dge.pval <- fread(glue::glue("{dgePrefix}/{contrast}_dge_genes.tsv"))


# ----------------- Save data -----------------
rds <- list(
  ds = an.intron,
  dge = dge.pval,
  dge_dsc = dge.dsc)

# make dir if not exists
if (!dir.exists(outPrefix)) {
  dir.create(outPrefix, recursive = TRUE)
}
out.f <- glue::glue("{outPrefix}/{contrast}_data.rds")
saveRDS(rds, out.f)

