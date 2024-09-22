# Coloc eQTL and sQTL using HyprColoc

if (interactive()) {
    args <- scan(
        text = "results/coloc/sqtl-eqtl/GTEx/Stomach/chr6.eqtl_sqtl_ids.txt
        results/coloc/sqtl-eqtl/GTEx/Stomach/chr6.eqtl_nominal.txt.gz
        results/coloc/sqtl-eqtl/GTEx/Stomach/chr6.sqtl_nominal.txt.gz
        ",
        what = character()
    )
    coloc_ids_f <- args[1]
    eqtl_nominal_f <- args[2]
    sqtl_nominal_f <- args[3]
} else {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) != 5) {
        stop("Usage: Rscript script.R coloc_ids_f eqtl_nominal_f sqtl_nominal_f out_txt out_rds")
    } else {
        coloc_ids_f <- args[1]
        eqtl_nominal_f <- args[2]
        sqtl_nominal_f <- args[3]
        out_txt <- args[4]
        out_rds <- args[5]
    }
}

suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
library(glue)
library(hyprcoloc)


#--- Load coloc eqtl-sqtl pids
print(glue("Loading coloc eqtl-sqtl pids from {coloc_ids_f}"))
coloc_ids <- fread(coloc_ids_f) %>% unique()
# coloc_ids[, id := glue("{pid}|{gid}", pid = phenotype_id, gid = gene_id)]
gene_ids <- coloc_ids$gene_id %>% unique()
names(gene_ids) <- gene_ids
print(glue("Number of genes for coloc: {length(gene_ids)}"))

#--- Load eQTL and sQTL nominal pass
summ_cols <- c(
               "pid", "pchr", "pstart", "pend", "pstrand",
               "nsnps", "distance", "snp_id", "snp_chr", "snp_start",
               "snp_end", "pval", "r2", "slope", "se",
               "top_snp"
)

eqtl_summ <- fread(eqtl_nominal_f, col.names = summ_cols)
sqtl_summ <- fread(sqtl_nominal_f, col.names = summ_cols)


#--- Useful functions

getMolQtlSumm <- function(nominal_f, coord) {
    # get eQtl or sQTL nominal pass (from qtltools) with tabix
    summ_cols <- c(
                   "pid", "pchr", "pstart", "pend", "pstrand",
                   "nsnps", "distance", "snp_id", "snp_chr", "snp_start",
                   "snp_end", "pval", "r2", "slope", "se",
                   "top_snp"
    )
    # coordinates must be 1-based, in format chr:start-end
    cmd <- glue("tabix {nominal_f} {coord}")
    summ <- fread(cmd = cmd, col.names = summ_cols)
    return(summ)
}

getTables <- function(gid, coloc_ids, eqtl, sqtl) {
    intron_ids <- coloc_ids[gene_id == gid, phenotype_id]
    print(glue("Prepare to coloc for {gid} and {paste(intron_ids, collapse = '\n')}\n"))

    # subset eqtl and sqtl summary stats, remove slope/se == 0
    eqtl <- eqtl[
        pid == gid & slope != 0 & se != 0,
        .(pid, snp_id, pval, slope, se)
    ]
    sqtls <- sqtl[
        pid %in% intron_ids & slope != 0 & se != 0,
        .(pid, snp_id, pval, slope, se)
    ] %>%
        split(by = "pid")

    # shared snps by the gene and all introns
    intron_shared_snps <- map(sqtls, ~ .x$snp_id) %>%
        purrr::reduce(intersect)
    if (length(intron_shared_snps) == 0) {
        print(glue("No shared snps for {gid}\n"))
        return(NULL)
    }
    shared_snps <- intersect(eqtl$snp_id, intron_shared_snps)
    N_shared_snps <- length(shared_snps)
    if (N_shared_snps == 0) {
        print(glue("No shared snps for {gid}\n"))
        return(NULL)
    } else {
        print(glue("Number of shared snps for {gid}: {N_shared_snps}\n"))
        
        # make sure the order of snps are the same before doing cbind
        eqtl <- eqtl[snp_id %in% shared_snps][order(snp_id)]
        sqtls <- map(sqtls, ~.x[snp_id %in% shared_snps] %>% .[order(snp_id)])
        mx.cols <- c(gid, map_chr(sqtls, ~.x$pid[[1]]))
        rowids = eqtl$snp_id

        betas <- cbind(eqtl$slope, map_dfc(sqtls, ~.x$slope))
        colnames(betas) <- mx.cols
        rownames(betas) <- rowids

        ses <- cbind(eqtl$se, map_dfc(sqtls, ~.x$se))
        colnames(ses) <- mx.cols
        rownames(ses) <- rowids

        return(list(betas = as.matrix(betas), ses = as.matrix(ses)))
    }
}


#--- Run coloc

coloc.results <- map(
    gene_ids,
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

# remove null results (no shared snps)
coloc.results <- compact(coloc.results)

# output coloc results
coloc.df <- imap_dfr(coloc.results, ~ pluck(.x, 1) %>% mutate(id = .y)) %>%
    filter(posterior_prob > 0)

if (nrow(coloc.df) > 0) {
    print(glue("Output coloc results to {out_txt} and {out_rds}"))
    fwrite(coloc.df, out_txt, sep = "\t", col.names = TRUE, quote = FALSE)
    write_rds(coloc.results, out_rds)
} else {
    print("No coloc found.")
}


