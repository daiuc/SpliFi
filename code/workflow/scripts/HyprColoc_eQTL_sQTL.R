# Coloc eQTL and sQTL using HyprColoc

if (interactive()) {
    args <- scan(
        text = "results/coloc/sqtl-eqtl/GTEx/Liver/chr1.eqtl_sqtl_ids.txt 
        results/coloc/sqtl-eqtl/GTEx/Liver/chr1.eqtl_nominal.txt.gz 
        results/coloc/sqtl-eqtl/GTEx/Liver/chr1.sqtl_nominal.txt.gz 
        ",
        what = character()
    )
    coloc_ids_f <- args[1]
    eqtl_nominal_f <- args[2]
    sqtl_nominal_f <- args[3]
} else {
    args <- commandArgs(trailingOnly = TRUE)
    coloc_ids_f <- args[1]
    eqtl_nominal_f <- args[2]
    sqtl_nominal_f <- args[3]
    out_f <- args[4]
    if (length(args) == 3) {
        stop("Rscript HyprColoc_eQTL_sQTL.R <coloc_ids> <eqtl_nominal> <sqtl_nominal> <output_file>")
    }
}

suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
library(glue)
library(hyprcoloc)


#--- Load coloc eqtl-sqtl pids
print(glue("Loading coloc eqtl-sqtl pids from {coloc_ids_f}"))
coloc_ids <- fread(coloc_ids_f) %>% unique
coloc_ids[, id := glue("{pid}|{gid}", pid = phenotype_id, gid = gene_id)]
print(glue("Number of coloc eqtl-sqtl pairs: {nrow(coloc_ids)}"))

#--- Load eqtl and sqtl summary
summ_cols <- c("pid", "pchr", "pstart", "pend", "pstrand",
               "nsnps", "distance", "snp_id", "snp_chr", "snp_start",
               "snp_end", "pval", "r2", "slope", "se",
               "top_snp")
eqtl_summ <- fread(eqtl_nominal_f, col.names = summ_cols)
sqtl_summ <- fread(sqtl_nominal_f, col.names = summ_cols)


#--- Useful functions

getTables <- function(coloc_id, coloc_ids, eqtl, sqtl) {

    gene_id <- coloc_ids[id == coloc_id, gene_id]
    intron_id <- coloc_ids[id == coloc_id, phenotype_id]

    eqtl <- eqtl[pid == gene_id & slope != 0 & se != 0,
                 .(pid, snp_id, pval, slope, se)]
    sqtl <- sqtl[pid == intron_id & slope != 0 & se != 0,
                 .(pid, snp_id, pval, slope, se)]
    shared_snps <- intersect(eqtl$snp_id, sqtl$snp_id)
    N_shared_snps <- length(shared_snps)
    if (N_shared_snps == 0) {
        print(glue("No shared snps for {coloc_id}"))
        return(NULL)
    } else {
        eqtl <- eqtl[snp_id %in% shared_snps]
        sqtl <- sqtl[snp_id %in% shared_snps]

        betas <- inner_join(
                    x = eqtl[, .(snp_id, slope)],
                    y = sqtl[, .(snp_id, slope)],
                    by = "snp_id",
                    suffix = c("_e", "_s")
        ) %>% arrange(snp_id) %>% column_to_rownames("snp_id")

        ses <- inner_join(
                    x = eqtl[, .(snp_id, se)],
                    y = sqtl[, .(snp_id, se)],
                    by = "snp_id",
                    suffix = c("_e", "_s")
        ) %>% arrange(snp_id) %>% column_to_rownames("snp_id")
        names(betas) <- c("eqtl", "sqtl")
        names(ses) <- c("eqtl", "sqtl")

        return(list(betas = as.matrix(betas), ses = as.matrix(ses)))
    }
}



#--- Run coloc

coloc.results <- map(
    coloc_ids$id,
    \(id) {
        mx <- getTables(id, coloc_ids, eqtl_summ, sqtl_summ)
        if (!is.null(mx)) {
            trait.names <- colnames(mx$betas)
            snp.ids <- rownames(mx$betas)
            res <- hyprcoloc(mx$betas, mx$ses,
                    trait.names = trait.names,
                    snp.id = snp.ids,
                    snpscores = T,
                    bb.selection = "alignment"
                    )
            return(res)
        } else {
            return(NULL)
        }
    }
)
names(coloc.results) <- coloc_ids$id

# remove null results (no shared snps)
coloc.results <- compact(coloc.results)

# output coloc results
coloc.df <- imap_dfr(coloc.results, ~pluck(.x, 1) %>% mutate(id = .y)) %>% 
    filter(posterior_prob > 0)

if (nrow(coloc.df) > 0) {
    print(glue("Output coloc results: {out_f}"))
    fwrite(coloc.df, out_f, sep = "\t", col.names = TRUE, quote = FALSE)
} else {
    print("No coloc found.")
}


