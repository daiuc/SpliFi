if (interactive()) {
    args <- scan(
        text = "resources/GWAS/hg38/AD_hitloci_1en5.bed 
        resources/GWAS/hg38/AD.tsv.gz 
        /project/yangili1/cdai/annotations/hg38/gencode.v26.GRCh38.genes.csv 
        results/coloc/sqtl-eqtl/GTEx/Liver 
        ",
        what = character()
    )
    gwas_loci_f <- args[1]
    gwas_stats_f <- args[2]
    genes_f <- args[3]
    sqtl_eqtl_dir <- args[4]
} else {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) != 5) {
        stop("Rscript HyprColoc_gwas_eQTL_sQTL.R <gwas_loci> <gwas_stats> <genes> <sqtl_eqtl_dir> <output_dir>")
    }
    gwas_loci_f <- args[1]
    gwas_stats_f <- args[2]
    genes_f <- args[3]
    sqtl_eqtl_dir <- args[4]
    out_dir <- args[5]
}


suppressMessages(library(tidyverse))
suppressMessages(library(glue))
suppressMessages(library(data.table))
suppressMessages(library(GenomicRanges))
suppressMessages(library(hyprcoloc))



#--- Load GWAS hit Loci
print(glue("Loading GWAS loci from {gwas_loci_f}"))
loci <- rtracklayer::import.bed(gwas_loci_f)
mcols(loci) <- NULL
gtrait = basename(gwas_loci_f) %>% str_split('_', simplify = TRUE) %>% .[1]


#--- with GWAS loci, get eQTLs and sQTLs (pids) for coloc analysis
print(glue("Loading sQTL-eQTL coloc pairs from {sqtl_eqtl_dir}"))
tissue = str_split(sqtl_eqtl_dir, '/', simplify = TRUE) %>% .[, 5]
pheno_ids = dir(glue('../code/results/coloc/sqtl-eqtl/GTEx/{tissue}'), '.*eqtl.*sqtl.*id.*txt', full.names = TRUE)
pheno_ids = dir(glue("{sqtl_eqtl_dir}"),'.*eqtl.*sqtl.*id.*txt', full.names = TRUE) %>% 
    naturalsort::naturalsort()
pheno_ids = map_dfr(pheno_ids, ~fread(.x)) %>% unique
colnames(pheno_ids) = c('intron_id', 'gene_id')

# Load gene annotations, used to get molQTL loci
genes = fread(genes_f)
genes = genes[feature == 'gene', .(seqname, start, end, strand, gene_id)] %>% unique()
genes[, gene_id := str_remove(gene_id, '\\..*')]
# only keep genes thave have molQTLs
genes = genes[gene_id %in% pheno_ids$gene_id]

# intersect genes with hit loci, keep only ones overlaps with hit loci
genes = makeGRangesFromDataFrame(genes, keep.extra.columns = TRUE)
# first remove genes that do not have any overlap with hit loci at all
genes.o = subsetByOverlaps(genes, loci, minoverlap = 10) 
# then use dynamic minoverlap
# require at least 10% of gene length overlaps with hit loci
minOlaps = round(width(genes.o) * 0.1)
genes.o = split(genes.o, seq_along(genes.o)) %>% 
    map2(minOlaps, ~subsetByOverlaps(.x, loci, minoverlap = .y))
# keep only ones that have overlaps
genes.o = genes.o[map_lgl(genes.o, ~length(.x) > 0)]
# unlist to Granges
names(genes.o) = NULL
genes.o = do.call(c, genes.o)
genes.o = sort(genes.o) # genes that overlap with hit loci

# these are the genes/introns that coloc with hit GWAS loci
# use these genes to select eqtl and sqtl
coloc.gene_ids = intersect(genes.o$gene_id, pheno_ids$gene_id)
coloc.intron_ids = pheno_ids[gene_id %in% coloc.gene_ids]$intron_id %>% unique

# intron id and gene id used for coloc, used to subset sqtl, eqtl
# lookup table: cols are intron_id, gene_id
coloc_phenos = pheno_ids[gene_id %in% coloc.gene_ids]

print(glue("Coloc analysis include: {gtrait}, {uniqueN(coloc_phenos$gene_id)} eQTLs and {uniqueN(coloc_phenos$intron_id)} sQTLs"))


#--- func for each locus
getSummaryByLocus = function(locus, gene_intron_lookup, gwas_stats_f, sqtl_eqtl_dir) {
    # given a locus (GRange),
    # get gwas stats for all snps within locus
    # get eqtl stats for all snps within locus, separate by gene
    # get sqtl stats for all snps within locus, separate by gene and intron

    print(glue("Getting summary stats for locus {as.character(locus)}"))
    chrom = as.character(seqnames(locus))
    gwas_stats_f = gwas_stats_f
    gtrait = basename(gwas_stats_f) %>% str_split('\\.', simplify = TRUE) %>% .[1]
    eqtl_stats_f <- glue("{sqtl_eqtl_dir}/{chrom}.eqtl_nominal.txt.gz")
    sqtl_stats_f <- glue("{sqtl_eqtl_dir}/{chrom}.sqtl_nominal.txt.gz")

    gwas.cols = fread(cmd = glue("zcat {gwas_stats_f} | head -n 2")) %>% colnames
    qtl.cols = c('pid', 'pchr', 'pstart', 'pend',
                 'pstrand', 'nsnps', 'distance', 'snp_id',
                 'snp_chr', 'snp_start', 'snp_end', 'P',
                 'r2', 'BETA', 'SE', 'top_snp')

    gwas_stats = suppressWarnings(fread(cmd = glue("tabix {gwas_stats_f} {as.character(locus)}"), col.names = gwas.cols))
    eqtl_stats = suppressWarnings(fread(cmd = glue("tabix {eqtl_stats_f} {as.character(locus)}"), col.names = qtl.cols))
    sqtl_stats = suppressWarnings(fread(cmd = glue("tabix {sqtl_stats_f} {as.character(locus)}"), col.names = qtl.cols))

    if (nrow(gwas_stats) == 0) {
        print(glue("No GWAS stats for locus {as.character(locus)}, exiting."))
        return(NULL)
    } else if (nrow(eqtl_stats) == 0) {
        print(glue("No eQTL stats for locus {as.character(locus)}, exiting."))
        return(NULL)
    } else if (nrow(sqtl_stats) == 0) {
        print(glue("No sQTL stats for locus {as.character(locus)}, exiting."))
        return(NULL)
    }

    # consistent snp_id to join on
    gwas_stats[, snp_id := paste0(`#CHR`, ':', BP)]
    eqtl_stats[, snp_id := paste0(snp_chr, ':', snp_start)]
    sqtl_stats[, snp_id := paste0(snp_chr, ':', snp_start)]

    # sometimes gwas has duplicates of snp_id
    gwas_stats[, rk := rank(P, ties.method = 'first'), by = snp_id]
    gwas_stats = gwas_stats[rk == 1]

    genes = base::intersect(eqtl_stats$pid, gene_intron_lookup$gene_id)
    print(glue("Locus {as.character(locus)}, width: {width(locus)}, overlaps 0 genes, exiting."))
    if (length(genes) == 0) return(NULL)
    names(genes) = genes
    print(glue("Locus {as.character(locus)}, width: {width(locus)}, overlaps {length(genes)} genes: {paste0(genes, collapse = ', ')}"))

    # for each gene, find shared SNPs between gwsa, eqtl, sqtl
    snps = map(
               genes,
               \(gid) {
                   # first shared snps across sqtls of the same gene
                   introns = base::intersect(
                                gene_intron_lookup[gene_id == gid, intron_id],
                                sqtl_stats$pid)
                   snps.sqtl = map(introns, ~sqtl_stats[pid == .x, snp_id]) %>% compact(.)
                   if (length(snps.sqtl) == 0) return(NULL)
                   snps.sqtl = purrr::reduce(snps.sqtl, base::intersect)
                   if (length(snps.sqtl) < 100) return(NULL)

                   # then shared snps across eqtl and sqtls of the gene, and gwas
                   snps.eqtl = eqtl_stats[pid == gid, snp_id] %>% unique
                   snps.gwas = gwas_stats$snp_id %>% unique
                   snps = list(snps.eqtl, snps.sqtl, snps.gwas) %>%
                       purrr::reduce(., base::intersect)
                   if (length(snps) >= 100) return(snps)
               }
    )
    snps = compact(snps) # list of shared snps for each gene
    if (length(snps) == 0) return(NULL)

    genes = names(snps) # genes with 100+ shared snps
    snps = map(snps, sort)
    print(glue("With {length(snps)} genes having 100+ shared snps across GWAS, eQTL, sQTL."))
    print(map_int(snps, length))

    # get gwas, eqtl, sqtl stats for shared snps
    data = imap(snps, \(shared_snps, genei) {

        shared_snps = shared_snps[order(shared_snps)]
        print(glue("have {length(shared_snps)} shared snps for {genei}"))
        introns = base::intersect(
                     gene_intron_lookup[gene_id == genei, intron_id],
                     sqtl_stats$pid)
        names(introns) = introns
        data.sqtl = map(introns,
                        \(intron) {
                            sqtl_stats[pid == intron & snp_id %in% shared_snps
                                ][, .(BETA, SE, P, rk = rank(P, ties.method = 'first')), by = snp_id
                                ][rk == 1, .(snp_id, BETA, SE)]
                        }
                        ) %>%
                    compact
        data.sqtl = map(data.sqtl, ~.x[order(snp_id)]) # ensure same order

        data.eqtl = eqtl_stats[pid == genei & snp_id %in% shared_snps
                        ][, .(snp_id, BETA, SE, P, rk = rank(P, ties.method = 'first')), by = snp_id
                        ][rk == 1, .(snp_id, BETA, SE)
                        ][order(snp_id)] # ensure smame order and unique snp_id

        data.gwas = gwas_stats[snp_id %in% shared_snps, .(snp_id, BETA, SE)][order(snp_id)]

        # beta.mx = cbind(data.gwas$BETA,
        #                 data.eqtl$BETA,
        #                 do.call(cbind, map(data.sqtl, ~.x$BETA))
        #                 )
        
        beta.mx = do.call("c", 
                          list(list(data.gwas$BETA),
                               list(data.eqtl$BETA),
                               map(data.sqtl, ~.x$BETA))
                          ) %>% as.data.frame
        colnames(beta.mx) = c(gtrait, genei, names(data.sqtl))
        rownames(beta.mx) = shared_snps

        se.mx = do.call("c", 
                        list(list(data.gwas$SE),
                             list(data.eqtl$SE),
                             map(data.sqtl, ~.x$SE))
                        ) %>% as.data.frame
        colnames(se.mx) = c(gtrait, genei, names(data.sqtl))
        rownames(se.mx) = shared_snps

        return(list(beta = as.matrix(beta.mx), se = as.matrix(se.mx)))
    })

    return(data)
}


runHyprColoc = function(x) {
    # x: list of 2 matrcis, beta and se
    # run hyprcoloc on each dataframe
    beta = x$beta
    se = x$se

    # remove rows with any zeros because hyprcoloc does not allow zeros
    non_zeros1 = rownames(beta)[rowSums(beta != 0) == ncol(beta)]
    non_zeros2 = rownames(se)[rowSums(se != 0) == ncol(se)]
    keep_rows = intersect(non_zeros1, non_zeros2)
    beta = beta[keep_rows, ] %>% as.matrix
    se = se[keep_rows, ] %>% as.matrix

    trait.names = colnames(beta)
    snp.ids = rownames(beta)
    res = hyprcoloc(beta, se,
                    trait.names = trait.names,
                    snp.id = snp.ids,
                    snpscores = TRUE,
                    bb.selection = "alignment"
                    )

    return(res)
}

resToDf = function(res, locus, gtrait) {
    # res: hyprcoloc result of a locus, a list
    # return a dataframe

    df = imap_dfr(res, ~pluck(.x, "results") %>% 
                  as.data.frame %>%
                  mutate(coloc_id = glue("{gtrait}|{as.character(locus)}|{.y}")
                         )
                  ) %>%
        filter(posterior_prob > 0)
    return(df)
}


#--- Run hyprcoloc for each locus

loci.l = split(loci, seq_along(loci))
names(loci.l) = NULL
coloc_results = map(
    loci.l,
    \(locus) {
        inputdata = getSummaryByLocus(locus, pheno_ids, gwas_stats_f, sqtl_eqtl_dir)
        if (is.null(inputdata)) return(NULL)
        res = map(inputdata, runHyprColoc)
        res.df = resToDf(res, locus, gtrait) # extract results to df
        if (nrow(res.df) > 0) {
            return(list(rds = res, df = res.df))
        }
    }
)
names(coloc_results) = as.character(loci)
coloc_results = compact(coloc_results)


#--- write coloc results in txt and rds

if (length(coloc_results) == 0) {
    stop("Found 0 loci with coloc results, exiting.")
}

print(glue("Found {length(coloc_results)} loci with coloc results"))
print(glue("Writing coloc results to {out_dir}/{gtrait}-eqtl-sqtl.coloc.tsv and {out_dir}/{gtrait}-eqtl-sqtl.coloc.rds"))

map_dfr(coloc_results, ~pluck(.x, 'df'))  %>% 
    fwrite(glue("{out_dir}/{gtrait}-eqtl-sqtl.coloc.tsv"), sep = "\t", col.names = TRUE, quote = FALSE)
map(coloc_results, ~pluck(.x, 'rds')) %>% 
    write_rds(glue("{out_dir}/{gtrait}-eqtl-sqtl.coloc.rds"), compress = 'gz')
print("Done.")

