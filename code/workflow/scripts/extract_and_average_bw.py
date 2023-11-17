import argparse
import numpy as np
import pandas as pd
import pysam
import pyBigWig
import glob
from jinja2 import Template
from os import path

def get_samples(sampleFile):
    smp = []
    with open(sampleFile)  as inF:
        for l in inF:
            smp.append(l.strip())
    return(smp)

def split_snpName(snpName):
    chrm, snpCoord = snpName.split("-")
    snpCoord = int(snpCoord)
    return(chrm, snpCoord)


def get_genotype(VCFFileName, chrm, snpCoord):
    vcfFile = pysam.VariantFile(VCFFileName)
    allSnp = vcfFile.fetch(chrm, snpCoord - 1, snpCoord)
    for snp in allSnp:
        if snp.pos == snpCoord:
            return(snp)


def get_bw_dict(genotype, bwFilePattern, samples):
    chrm = genotype.contig
    snpCoord = genotype.pos
    bwDict = {0: [], 1: [], 2:[]}
    smpDict = {0: [], 1: [], 2:[]}
    for smp in genotype.samples.keys():
        if samples is not None and smp in samples or samples is None:
            tmpBwFile = bwFilePattern.format(smp)
            tmpGeno = genotype.samples[smp]['GT']
            if tmpGeno == (0,0):
                bwDict[0].append(tmpBwFile)
                smpDict[0].append(smp)
            elif tmpGeno in [(0, 1), (1, 0)]:
                bwDict[1].append(tmpBwFile)
                smpDict[1].append(smp)
            elif tmpGeno == (1, 1):
                bwDict[2].append(tmpBwFile)
                smpDict[2].append(smp)
    return((bwDict, smpDict))


def get_bw_dict_from_file(fileName):
    bwDict = {}
    with open(fileName, "rt") as tmpFile:
        for ln in tmpFile:
            lnVec = ln.rstrip().split()
            try:
                bwDict[lnVec[1]].append(lnVec[0])
            except KeyError:
                bwDict[lnVec[1]] = [lnVec[0]]
    return(bwDict)


def args_parser():
    parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--VCF",
            help= "gzipped and tbi indexed vcf file with sample genotypes",
            metavar="<filename>",
            default=None)
    parser.add_argument("--SNP",
            help= "SNP position like chr1:12345. Position should be one-based, same as in vcf",
            metavar="<CHR>:<POS>",
            default=None)
    parser.add_argument("--averageBy",
            help= """If not using a VCF file, supply this with a tab-delimited file without header that has bw file name as the first column,
            and grouping variable as the second column.
            bw files will be averaged by the grouping variable.""",
            metavar="<filename>")
    parser.add_argument("--sampleFile",
            help= """
            Optional. A file containing each sample per line for teh samples to use.
            This should be a subset of all samples in the VCF file and BigWig files.
            If not used, then all samples are considered.""",
            metavar="<filename>",
            default=None)
    parser.add_argument("--Chr",
            help= "Chromosome name if SNP is not specified.",
            metavar="<CHR>",
            default=None)
    parser.add_argument("--regionStart",
            help= "Start position for region",
            metavar="<integer>")
    parser.add_argument("--regionEnd",
            help= "End position for region",
            metavar="<integer>")
    parser.add_argument("--extend",
            help= "Extend the cluster on both ends to help plotting",
            metavar="<integer>")
    parser.add_argument("--bwFilePattern",
            help= """
            A string in the form of data/bamCoverage/\{\}.bw where {} should be the place where sample name is.
            Sample name must be present in the VCF file.""",
            metavar="<string>")
    parser.add_argument("--template",
            help= """
            A jinja template for making pyGenomeTracks .ini file.""",
            default=None)
    parser.add_argument("--OutputPrefix",
            help="Prefix for all output files (default: %(default)s)",
            metavar="<output prefix>",
            default="./")
    return(parser.parse_args())


args = args_parser()
if args.sampleFile is not None:
    smps = get_samples(args.sampleFile)
else:
    smps = None

if args.SNP is not None:
    chrm, snpCoord = split_snpName(args.SNP)
    snpBedFile = args.OutputPrefix + f"{args.SNP}.bed"
    with open(snpBedFile, "wt") as vline:
        vline.write(f"{chrm}\t{snpCoord - 1}\t{snpCoord}")

if args.VCF is not None:
    geno = get_genotype(args.VCF, chrm, snpCoord)
    bwFileDict, smpDict = get_bw_dict(geno, args.bwFilePattern, smps)
    # Write smp table
    smpGrpName = args.OutputPrefix + f"smp_{args.SNP}.txt"
    with open(smpGrpName, "wt") as smpGrp:
        for allele in smpDict.keys():
            for smp in smpDict[allele]:
                smpGrp.write(f"{allele}\t{smp}\n")
elif args.averageBy is not None:
    bwFileDict = get_bw_dict_from_file(args.averageBy)

if args.Chr is not None:
    chrm = args.Chr

regionStart = int(args.regionStart)
regionEnd = int(args.regionEnd)
extend = int(args.extend)

avgBwDict = {}
# maxSmpNumber = max(len(bwFileDict.values()))

for allele, bwFiles in bwFileDict.items():
    avgBwDict[allele] = []
    for bwFile in bwFiles:
        if not path.isfile(bwFile):
            continue
        with pyBigWig.open(bwFile, "rt") as bwIn:
            bwHeader = list(bwIn.chroms().items())
            try:
                tmpValue = bwIn.values(chrm, regionStart - extend, regionEnd + extend, numpy=True)
            except RuntimeError:
                tmpValue = bwIn.values("chr" + chrm, regionStart - extend, regionEnd + extend, numpy=True)
            tmpValue = np.nan_to_num(tmpValue, nan=0.0)
            avgBwDict[allele].append(tmpValue)

# Write average bigwig files
avgFileDict = {} # used for jinja template
for allele, values in avgBwDict.items():
    if values == []:
        continue
    outName = args.OutputPrefix + f"avg_{args.SNP}_{regionStart}_{regionEnd}_{allele}.bw"
    avgFileDict[allele] = outName
    with pyBigWig.open(outName, "wt") as bwAvg:
        valueOut = np.mean(np.array(values), axis=0)
        bwAvg.addHeader(bwHeader)
        try:
            bwAvg.addEntries(chrm, regionStart - extend, values=valueOut, span=1, step=1)
        except RuntimeError:
            bwAvg.addEntries("chr"+chrm, regionStart - extend, values=valueOut, span=1, step=1)

if args.template != None:
    colorDict = dict(zip(geno.alleles, ["red", "blue"]))
    with open(args.template, "rt") as tplt:
        template_file_contents = tplt.read()
        template = Template(template_file_contents)
        fill_template = template.render(
            yMax=yMax,
            alleles=geno.alleles,
            bwFileDict=bwFileDict,
            avgFileDict=avgFileDict,
            colorDict=colorDict,
            nBins = args.nBins,
            snpBed = snpBedFile
        )

    with open(args.OutputPrefix + f"{chrm}_{snpCoord}.ini", "wt") as iniOut:
        iniOut.write(fill_template + "\n")
