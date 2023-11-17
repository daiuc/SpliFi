#!/usr/bin/env Rscript

# written by bjf79

#Use hard coded arguments in interactive R session, else use command line args
if (interactive()){
  args <- scan(text=
                 "/project2/yangili1/cdai/mRNA-editing/YRI/Results/hs38/edQTL_FixedPCs/cis_10000/permutations_pass_all_chr.txt.gz scratch/Qvals.txt.gz   ", what='character')
} else{
  args <- commandArgs(trailingOnly=TRUE)
}

FileIn <- args[1]
FileOut <- args[2]

suppressMessages(library(data.table))
suppressMessages(library(qvalue))


data.in = fread(FileIn)

if (ncol(data.in) == 20) { # QTLtools 1.3 output
  cnames <- c(
    "phenotype_id", "phenotype_chr", "phenotype_start", "phenotype_end",
    "phenotype_strand", "num_variants", "best_nom_dist", "best_genotype_id",
    "best_genotype_chr", "best_genotype_start", "best_genotype_end", "dof_true",
    "dof_est", "beta_ml1", "beta_ml2", "pval_nom",
    "pval_r2", "slope", "pval_emp", "pval_adj")
  names(data.in) <- cnames
} else if (ncol(data.in) == 19) { # QTLtools 1.2 output
  cnames <- c(
    "phenotype_id", "phenotype_chr", "phenotype_start", "phenotype_end",
    "phenotype_strand", "num_variants", "best_nom_dist", "best_genotype_id",
    "best_genotype_chr", "best_genotype_start", "best_genotype_end", "dof_true",
    "dof_est", "beta_ml1", "beta_ml2", "pval_nom",
    "slope", "pval_emp", "pval_adj")
  names(data.in) <- cnames
}

data.in$q <- signif(qvalue(data.in$pval_adj)$qvalues, 5)
fwrite(na.omit(data.in), FileOut, sep = " ")
