# Process QTLtools nominal pass results

if (interactive()) {
    args <- scan(
        text = "results/coloc/sqtl-eqtl/GTEx/Liver/chr1.eqtl_nominal.txt ",
        what = character()
    )
    nominal_f <- args[1]
} else {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) != 2) {
        stop("Rscript FormatQtltoolsNominalResult.R <nominal_file> <output_file>")
    }
    nominal_f <- args[1]
    out_f <- args[2]
}


suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
library(glue)

#--- load nominal pass qtls
print(glue("Load nominal pass qtls: \n{nominal_f}"))
nominal <- fread(nominal_f, header = FALSE)

# col1 pid
# col2 pchr
# col3 pstart
# col4 pend
# col5 pstrand
# col6 number of snps
# col7 distance
# col8 snp id
# col9 snp chr
# col10 snp start
# col11 snp end
# col12 nominal pval
# col13 r2
# col14 slope
# col15 slope SE
# col16 1/0 to indicate top snp

if (ncol(nominal) != 16) {
    stop(glue("The input file {nominal_f} does not have 16 fields"))
}

if (str_detect(nominal$V9[[1]], "chr")) {
    chroms <- c(glue("chr{1:22}"), "chrX", "chrY")
    nominal[, V9 := factor(V9, levels = chroms)]
}

nominal <- nominal[order(V9, V10, V11)]

#--- output
print(glue("Output tabix compatible: \n{out_f}"))
fwrite(nominal, out_f, sep = "\t", col.names = FALSE, quote = FALSE)
