#!/usr/bin/env Rscript

# written by bjf79

# Use hard coded arguments in interactive R session, else use command line args
if (interactive()) {
    args <- scan(
        text =
            "/project/yangili1/cdai/SpliFi/code/results/eqtl/GTEx/Testis/perm/chr18.txt scratch/Qvals.txt.gz   ", what = "character"
    )
} else {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) < 2) {
        stop("Rscript AddQvalueToQTLtoolsOutput.R <input file> <output file>")
    }
}

FileIn <- args[1]
FileOut <- args[2]

suppressMessages(library(data.table))
suppressMessages(library(qvalue))


data.in <- fread(FileIn)

if (ncol(data.in) == 20) { # QTLtools 1.3 output
    cnames <- c(
        "phenotype_id", "phenotype_chr", "phenotype_start", "phenotype_end",
        "phenotype_strand", "num_variants", "best_nom_dist", "best_genotype_id",
        "best_genotype_chr", "best_genotype_start", "best_genotype_end", "dof_true",
        "dof_est", "beta_ml1", "beta_ml2", "pval_nom",
        "pval_r2", "slope", "pval_emp", "pval_adj"
    )
    names(data.in) <- cnames
} else if (ncol(data.in) == 19) { # QTLtools 1.2 output
    cnames <- c(
        "phenotype_id", "phenotype_chr", "phenotype_start", "phenotype_end",
        "phenotype_strand", "num_variants", "best_nom_dist", "best_genotype_id",
        "best_genotype_chr", "best_genotype_start", "best_genotype_end", "dof_true",
        "dof_est", "beta_ml1", "beta_ml2", "pval_nom",
        "slope", "pval_emp", "pval_adj"
    )
    names(data.in) <- cnames
}

q <- try(qvalue(data.in$pval_adj), silent = TRUE)
if (inherits(q, "try-error")) {
    message("Setting qvalue parameter lambda=0 and rerun qvalue.")
    q <- qvalue(data.in$pval_adj, lambda = 0)
}
data.in$q <- signif(q$qvalues, 5)
fwrite(na.omit(data.in), FileOut, sep = " ")
