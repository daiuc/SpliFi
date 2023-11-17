### Clean up phenotype data from leafcutter prep pheno output


library(tidyverse)
library(data.table)


if (!interactive()) {
  print("Running in script mode")
  myargs = commandArgs(trailingOnly = T)

  if (length(myargs) < 2) {
    print("Error. Enter command 'Rscript path/to/Make_pheno_bed.R pheno_file outfile'")
    quit(save = "no", status = 1)
  } else {
    print(myargs)
    input_file = myargs[[1]]
    out_file = myargs[[2]]
  }

} else {
  print("Running in interactive/dev mode")

  input_file = "results/Geuvadis/noise_pheno/Geuvadis_pheno_extract.count.gz.phen_chr22"
  sample_pat = "[A-Z]{2}[0-9]{5,7}" # sample name regex pattern
}

# pheno type table
pheno = fread(input_file)

sample_pat = "[A-Z]{2}[0-9]{5,7}" # sample name regex pattern
datacols = names(pheno)[str_detect(names(pheno), sample_pat)] %>%
  str_extract(sample_pat)

idcols = names(pheno)[!str_detect(names(pheno), sample_pat)]

names(pheno) = c(idcols, datacols)

pheno[, `:=`(
  pid = str_extract(ID, "chr[0-9]+\\:[0-9]+\\:[0-9]+\\:clu_[0-9]+"),
  gid = str_extract(ID, "chr[0-9]+\\:[0-9]+\\:[0-9]+\\:clu_[0-9]+"),
  strand = str_extract(ID, "[\\+\\-]")
)]


selectcols = c(c("#Chr", "start", "end", "pid", "gid", "strand"), datacols)
pheno = pheno[order(`#Chr`, start, end), ..selectcols]

fwrite(pheno, out_file, sep = "\t")
print("finished writing to file.")
