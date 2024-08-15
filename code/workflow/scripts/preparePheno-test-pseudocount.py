#!/usr/bin/env python
"""
Modified from leafcutter 1 preparation script
"""


import gzip
import sys
from optparse import OptionParser

import numpy as np
from scipy.stats import norm, rankdata
from sklearn import preprocessing
from sklearn.decomposition import PCA

# import scipy as sc


def qqnorm(x):
    n = len(x)
    a = 3.0 / 8.0 if n <= 10 else 0.5
    return norm.ppf((rankdata(x) - a) / (n + 1.0 - 2.0 * a))


def stream_table(f, ss=" "):
    fc = "chrom"  # first column, different from leafcutter1
    while fc == "chrom":
        fc = f.readline().strip()
        head = fc.split(ss)

    for ln in f:
        ln = ln.strip().split(ss)
        attr = {}

        for i in range(len(head)):
            try:
                attr[head[i]] = ln[i]
            except:
                break
        yield attr


def get_chromosomes(ratio_file):
    """Get chromosomes from table. Returns set of chromosome names"""
    try:
        open(ratio_file)
    except:
        sys.stderr.write(f"Can't find {ratio_file}..exiting\n")
        return
    sys.stderr.write("Parsing chromosome names...\n")
    chromosomes = set()
    with gzip.open(ratio_file, "rt") as f:
        f.readline()
        for line in f:
            chromosomes.add(line.split(":")[0])
    return chromosomes


def get_blacklist_chromosomes(chromosome_blacklist_file):
    """
    Get list of chromosomes to ignore from a file with one blacklisted
    chromosome per line. Returns list. eg. ['X', 'Y', 'MT']
    """
    if chromosome_blacklist_file:
        with open(chromosome_blacklist_file, "r") as f:
            return f.read().splitlines()
    else:
        return ["X", "Y", "MT", "chrX", "chrY", "chrM"]


def main(
    ratio_file: str,
    sample_file: str,
    chroms: set,
    blacklist_chroms: list,
    pcs=50,
    outPrefix="out",
):

    fout = {}
    try:
        open(ratio_file)
    except:
        sys.stderr.write(f"Can't find {ratio_file}..exiting\n")
        return

    sys.stderr.write("Starting...\n")
    for i in chroms:
        fout[i] = open(outPrefix + ".phen_" + i, "w")  # phenotype out files
        fout_ave = open(outPrefix + ".ave", "w")
    valRows, valRowsnn, geneRows = [], [], []
    header = (
        gzip.open(ratio_file, "rt").readline().split()[1:]
    )  # sample names in counts file
    header_clean = [
        x.split(".")[0] for x in header
    ]  # remove file extension from sample names

    if sample_file:
        try:
            sample_names = [x.strip() for x in open(sample_file, "r").readlines()]
            sys.stderr.write(
                f"Filtering using {len(sample_names)} sample names from {sample_file}...\n"
            )
            sys.stderr.write(f"Counts file has {len(header)} samples, ")
            old_header = header
            # remove samples not in vcf sample list
            header = [x for x in header if x.split(".")[0] in sample_names]
            header_clean = [
                x.split(".")[0] for x in header
            ]  # remove file extension from sample names
            if len(header) == 0:

                sys.stderr.write(
                    f"Error: no sample names in {sample_file} match the sample names in the counts file. Exiting.\n"
                )
                exit(0)
            elif len(header) < len(old_header):
                lost_header = [x for x in old_header if x not in header]
                sys.stderr.write(
                    f"dropping {len(lost_header)} samples not in {sample_file}:\n{' '.join(lost_header)}\n"
                )

            sys.stderr.write(
                f"keeping {len(header)} samples matched in {sample_file}\n"
            )

        except:
            sys.stderr.write(f"Can't find {sample_file}..not filtering samples\n")
            sample_file = None

    for i in fout:  # write cleaned header to phenotype out files
        fout[i].write(
            "\t".join(["#Chr", "start", "end", "pid", "gid", "strand"] + header_clean)
            + "\n"
        )

    for dic in stream_table(gzip.open(ratio_file, "rt"), " "):

        if sample_file:
            dic = {
                k: v for k, v in dic.items() if k in ["chrom"] + header
            }  # header is intersect of samples in count files and vcf sample list

        chrom = dic["chrom"]
        chr_ = chrom.split(":")[0].replace("chr", "")
        if chr_ in blacklist_chroms:
            continue
        NA_indices, valRow, aveReads = [], [], []
        tmpvalRow = []

        i = 0
        for sample in header:
            try:
                count = dic[sample]
            except:
                print([chrom, len(dic)])
                print(dic)
            num, denom = count.split("/")
            if float(denom) < 1:
                count = "NA"
                tmpvalRow.append("NA")
                NA_indices.append(i)
            else:
                # addon = 0.5 # original pseudocount method
                addon = float(denom) * 0.01 if float(denom) < 100 else 0.01
                count = (float(num) + addon) / (float(denom) + addon)
                tmpvalRow.append(count)
                aveReads.append(count)

        # If ratio is missing for over 40% of the samples, skip
        if tmpvalRow.count("NA") > len(tmpvalRow) * 0.4:
            continue

        ave = np.mean(aveReads)

        # Set missing values as the mean of all values
        for c in tmpvalRow:
            if c == "NA":
                valRow.append(ave)
            else:
                valRow.append(c)

        # If there is too little variation, skip (there is a bug in fastqtl which doesn't handle cases with no variation)
        if np.std(valRow) < 0.005:
            continue

        chr_, s, e, clu, *_ = chrom.split(":")
        strand = clu.split("_")[2]
        gid = "."
        if len(valRow) > 0:
            fout[chr_].write(
                "\t".join([chr_, s, e, chrom, gid, strand] + [str(x) for x in valRow])
                + "\n"
            )
            fout_ave.write(
                " ".join(
                    [f"{chrom}"]
                    + [str(min(aveReads)), str(max(aveReads)), str(np.mean(aveReads))]
                )
                + "\n"
            )

            # scale normalize
            valRowsnn.append(valRow)
            valRow = preprocessing.scale(valRow)

            valRows.append(valRow)
            geneRows.append("\t".join([chr_, s, e, chrom, gid, strand]))
            if len(geneRows) % 5000 == 0:
                sys.stderr.write(f"Parsed {len(geneRows)} introns...\n")

    for i in fout:
        fout[i].close()

    # qqnorms on the columns
    matrix = np.array(valRows)
    for i in range(len(matrix[0, :])):
        matrix[:, i] = qqnorm(matrix[:, i])

    # write the qqnorm corrected tables
    fout = {}
    for i in chroms:
        fn = f"{outPrefix}.qqnorm_{i}"
        print("Outputting: " + fn)
        fout[i] = open(fn, "w")
        fout[i].write(
            "\t".join(["#Chr", "start", "end", "pid", "gid", "strand"] + header_clean)
            + "\n"
        )
    lst = []
    for i in range(len(matrix)):
        chrom, s = geneRows[i].split()[:2]

        lst.append(
            (
                chrom,
                int(s),
                "\t".join([geneRows[i]] + [str(x) for x in matrix[i]]) + "\n",
            )
        )

    lst.sort()
    for ln in lst:
        fout[ln[0]].write(ln[2])

    fout_run = open(f"{outPrefix}_prepare.sh", "w")

    for i in fout:
        fout[i].close()
        fout_run.write(f"bgzip -f {outPrefix.split('/')[-1]}.qqnorm_{i}\n")
        fout_run.write(f"tabix -p bed {outPrefix.split('/')[-1]}.qqnorm_{i}.gz\n")
    fout_run.close()

    sys.stdout.write(
        f"Use `sh {outPrefix}_prepare.sh' to create index for fastQTL (requires tabix and bgzip).\n"
    )

    # write sample names in to file
    fout_samples = open(f"{outPrefix}_names.txt", "w")
    sample_names = header_clean
    fout_samples.write("\n".join(sample_names))

    if pcs > 0:
        # matrix = np.transpose(matrix) # important bug fix (removed as of Jan 1 2018)
        pcs = min([len(header), pcs])
        pca = PCA(n_components=pcs)
        pca.fit(matrix)
        pca_fn = outPrefix + ".PCs"
        print("Outputting PCs: " + pca_fn)
        pcafile = open(pca_fn, "w")
        pcafile.write("\t".join(["id"] + header_clean) + "\n")
        pcacomp = list(pca.components_)

        for i in range(len(pcacomp)):
            pcafile.write("\t".join([str(i + 1)] + [str(x) for x in pcacomp[i]]) + "\n")

        pcafile.close()


if __name__ == "__main__":

    parser = OptionParser(usage="usage: %prog [-p num_PCs] input_perind.counts.gz")
    parser.add_option(
        "-p", "--pcs", dest="npcs", default=50, help="number of PCs output"
    )
    parser.add_option(
        "--ChromosomeBlackList",
        dest="cbl",
        default="",
        help="file of blacklisted chromosomes to exclude from analysis, one per line. If none is provided, will default to blacklisting X and Y",
    )
    parser.add_option(
        "--outPrefix",
        dest="outPrefix",
        default="leafcutter",
        help="prefix for output files",
    )
    parser.add_option(
        "--sampleFile",
        dest="sampleFile",
        default=None,
        help="List of sample names to filter counts file by. If none is provided, will use all samples in counts file",
    )
    (options, args) = parser.parse_args()
    if len(args) == 0:
        sys.stderr.write("Error: no ratio file provided...\n")
        exit(0)
    ratio_file, npcs, chrom_blacklist, outPrefix, sampleFile = (
        args[0],
        int(options.npcs),
        options.cbl,
        options.outPrefix,
        options.sampleFile,
    )
    chrom_list, chrom_blacklist = get_chromosomes(
        ratio_file
    ), get_blacklist_chromosomes(chrom_blacklist)
    chrom_list = chrom_list - set(chrom_blacklist)
    main(ratio_file, sampleFile, chrom_list, chrom_blacklist, npcs, outPrefix)
