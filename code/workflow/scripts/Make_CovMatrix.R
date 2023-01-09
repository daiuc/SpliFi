
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

if(interactive()) {
  print("Running in interactive mode")
  texts = "results/Geuvadis/noise_pheno/chr22_pheno.pca results/Genotype/Geuvadis/chr22_geno.pca 4"
  args = scan(text = texts, what = "character")
} else{
  print("Running in script mode")
  args <- commandArgs(trailingOnly=TRUE)
}

pheno_pca = args[[1]]
geno_pca = args[[2]]
geno_N_PCs = as.integer(args[[3]]) # number of pcs to choose from genotype pca
outfile = args[[4]]

pheno = fread(pheno_pca)
geno = fread(geno_pca)
geno = head(geno, geno_N_PCs)


cnames = names(pheno)
cnames[[1]] = "SampleID"
names(pheno) = cnames

pheno_samples = names(pheno)
pheno_samples = pheno_samples[2:length(pheno_samples)]

geno_samples = names(geno)
geno_samples = geno_samples[2:length(geno_samples)]



if (length(setdiff(geno_samples, pheno_samples)) == 0 &
    length(setdiff(pheno_samples, geno_samples)) == 0) {
      selectCols = c("SampleID", pheno_samples)
    } else {
      print("Genotype and Phenotype samples do not match! Exit.")
      quit(save = "no", status = 1)
    }

pheno = pheno[, ..selectCols]
geno = geno[, ..selectCols]

covMatrix = rbindlist(list(pheno, geno))

covMatrix = rename(covMatrix, id = SampleID)

fwrite(covMatrix, outfile, sep = "\t")