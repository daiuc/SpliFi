# most of these codes are from bjf79

suppressMessages(library(tidyverse))

if(interactive()) {
  print("Running in interactive mode")
  args = "results/Geuvadis/noise_pheno/chr22_pheno.bed.gz"
} else{
  print("Running in script mode")
  args <- commandArgs(trailingOnly=TRUE)
}

FileIn <- args[1]
FileOut <- args[2]

header = read_tsv(FileIn, n_max = 5)
ncols = ncol(header)

mat <- read_tsv(FileIn) %>%
    select(-c(1,2,3,5,6)) %>%
    column_to_rownames("pid") %>%
    as.matrix()


mat.permuted <- matrix(nrow=nrow(mat), ncol=ncol(mat))
for (i in 1:nrow(mat)){
  mat.permuted[i,] <- sample(mat[i,], size=ncol(mat))
}

pca.results <- prcomp(t(mat))
pca <- pca.results %>% summary() %>% pluck("importance") %>% t() %>% as.data.frame() %>%
  rownames_to_column("PC")

pca.permuted <- prcomp(t(mat.permuted)) %>% summary() %>% pluck("importance") %>% t() %>%
  as.data.frame() %>% rownames_to_column("PC")

merged <- bind_rows(list(pca=pca, pca.permuted=pca.permuted), .id="mat") %>%
  mutate(PC=as.numeric(str_replace(PC, "PC", "")))

#GetNumPCs
NumPCs <- merged %>%
  select(PC, Prop=`Proportion of Variance`, mat) %>%
  pivot_wider(names_from = "mat", values_from = "Prop") %>%
  filter(pca > 1.01*pca.permuted) %>% pull(PC) %>% max()

print("Num PCs for which variance explained is more than in permuted data")
print(NumPCs)

if (NumPCs > 15) {
  print("More than 20 PCs, set NumPC=15")
  NumPCs = 15
}

print("Writing out phenotype PCA results.")
pca.results$x[,1:NumPCs] %>% t() %>%
  round(5) %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  write_delim(FileOut)
