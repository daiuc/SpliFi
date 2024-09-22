"""
Prepare dataset required for plotting sashimi plots using pyGenomeTracks
"""

import glob
import os

import numpy as np
import pandas as pd
import pyBigWig as pw
from jinja2 import Template


def avgBigwig(bwfiles, grange):
    """Average multiple bigwig files in a specific region

    bwfiles : list of bigwig files (path).
    grange  : tuple. Genomic range, BED like 0-based coordinates,
              eg. ('chr1', 25101, 27101)

    return  : a dictionary of keys:
        - header, for pyBigWig to write as header
        - values, for pyBigwig to addentries as values
    """
    chrom, start, end = grange
    # compute much wider region than given range
    if end - start < 200000:
        start = max(5000, start - 200000)
        end = end + 200000
    values = []
    bwo = {}
    for bw in bwfiles:
        if not os.path.isfile(bw):
            continue
        with pw.open(bw, "rt") as b:
            header = list(b.chroms().items())
            vals = b.values(chrom, start, end, numpy=True)
            vals = np.nan_to_num(vals)
            values.append(vals)

    if values != [] and header != []:
        avgValues = np.mean(values, axis=0)
        bwo = {"header": header, "values": avgValues}
    return bwo


def writeBigwig(outFile, bigwig, grange):
    """Write a bigwig file

    outFile : str. Output bigwig file path
    bigwig  : dict. Dictionary with keys: header, values of a bigwig object
    grange  : tuple. Genomic range, BED like 0-based coordinates,
                eg. ('chr1', 25101, 27101)
    """
    chrom, start, end = grange
    if end - start < 200000:
        start = max(5000, start - 200000)

    with pw.open(outFile, "w") as b:
        b.addHeader(bigwig["header"])
        b.addEntries(chrom, start, values=bigwig["values"], span=1, step=1)
    print(f"Bigwig file written: {outFile}")


def readIntronTable(intronTableFile, sep=" ", SelectCols=[0]):
    """Read intron table file

    intronTableFile : str. Path to intron table file
                      eg. 'ds_perind_numers.counts.noise_by_intron.gz'

    return : pandas.DataFrame
    """
    introns = pd.read_csv(intronTableFile, sep=sep, usecols=SelectCols)
    introns = [
        (":".join(x[:-1]), f"{x[0]}:{x[3]}", x[4])
        for x in introns["chrom"].str.split(":")
    ]
    introns = {x[0]: x[2] for x in introns}

    return introns


def getSamplesByGroup(
    sampleGroupFile,
    groups,
    sep=" ",
    cols=["fname", "group"],
    header=None,
    splitBy=".",
    splitCapture=[1],
):
    """Read sample group file from diff splicing

    sampleGroupFile : str. Path to sample group file
                      eg. 'ds_sample_group.txt'
    sep             : str. Separator used in the file
    cols            : list. Column names to read from the file
    header          : int. Row number to use as header
    splitBy         : str. Split the file name by this character
    splitCapture    : list. Index to capture after splitting
    groups          : list. List of groups to select, eg. ['Liver', 'Kidney']


    return : dict. Dictionary of sample groups
    """
    samples = pd.read_csv(sampleGroupFile, sep=sep, header=header, names=cols)
    samples["sname"] = [x.split(splitBy)[splitCapture[0]] for x in samples["fname"]]

    group1, group2 = groups
    samples1 = samples[samples["group"] == group1]["sname"].tolist()
    samples2 = samples[samples["group"] == group2]["sname"].tolist()

    return {group1: samples1, group2: samples2}


def getBigWigFiles(bwDir, filePattern="*.bw"):
    """Get bigwig files from a directory

    bwDir      : str. Path to the directory containing bigwig files
    filePattern: str. File pattern to search for, eg: '*.bw'

    return     : dict. Dictionary of sample names and bigwig files
    """
    bwFiles = glob.glob(os.path.join(bwDir, filePattern))
    bwSamples = ["-".join(x.split("/")[-1].split("-")[:2]) for x in bwFiles]

    return {x[0]: x[1] for x in zip(bwSamples, bwFiles)}


def getPlottableBigWig(bw, ds):
    """Get plottable bigwig files
    bw : dict. key = sample name, value = bigwig file path
    ds : list. List of sample names in differential splicing.

    return: dict. Key = sample names in both bw and ds, value = bigwig file path
    """
    shared = set(bw.keys()).intersection(ds)
    return {x: bw[x] for x in shared}


def getPSI(effects, clu, introns, groups, minPSI=0.05):
    """Get PSI values from a file

    effects : dataframe.
    clu     : str. cluster id, eg. 'clu_261_+'
    introns: dict. Dictionary of introns
    groups  : list. List of groups, eg. ['Liver', 'Kidney']
    minPSI  : float. Minimum PSI value to output

    return : dict. Dictionary of links
    """
    effects["itype"] = [introns[x] for x in effects["intron"]]
    effects = effects[
        effects.intron.str.contains(clu)
    ].copy()  # use copy() to avoid SettingWithCopyWarning

    # make pyGenomeTracks compatible link tables
    effects.loc[:, "chrom1"] = effects.intron.str.split(":").str[0]
    effects.loc[:, "start1"] = effects.intron.str.split(":").str[1]
    effects.loc[:, "end1"] = effects["start1"]
    effects.loc[:, "chrom2"] = effects["chrom1"]
    effects.loc[:, "start2"] = effects.intron.str.split(":").str[2]
    effects.loc[:, "end2"] = effects["start2"]

    grp1, grp2 = groups
    outdfs = {
        grp1: effects[["chrom1", "start1", "end1", "chrom2", "start2", "end2", grp1]][
            effects[grp1] > minPSI
        ],
        grp2: effects[["chrom1", "start1", "end1", "chrom2", "start2", "end2", grp2]][
            effects[grp2] > minPSI
        ],
    }

    # links for UP junctions only
    upLinks = effects[effects.itype == "UP"].drop_duplicates()[
        ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]
    ]

    return {"sashimiLinks": outdfs, "upLinks": upLinks}


def getPlotRange(clu, introns, minRange=8000):
    """Get plot ranges for sashimi plot

    PlotClu : str. Cluster id, eg. 'clu_261_+'
    Introns : dict. Dictionary of introns

    return  : tuple. Genomic range, BED like 0-based coordinates,
              eg. ('chr1', 25101, 27101)
    """
    Keys = [x for x in introns.keys() if clu in x]
    KeysLs = [x.split(":") for x in Keys]
    CluChrom = KeysLs[0][0]
    CluStart = min([int(x[1]) for x in KeysLs])
    CluEnd = max([int(x[2]) for x in KeysLs])
    CluLen = CluEnd - CluStart
    if CluLen >= minRange:
        return (CluChrom, CluStart, CluEnd)
    elif CluLen < 5000:
        return (CluChrom, CluStart - 1500, CluEnd + 1500)
    else:
        d = (minRange - CluLen) // 2
        return (CluChrom, CluStart - d, CluEnd + d)


def writeIni(templateFile, outIniFile, trackFiles, nBins=700):
    """Generate ini file for plotting sashimi plots
    templateFile : str. Path to template ini file
    outIniFile   : str. Path of ini file to write into.
    trackFiles   : dict. Dictionary hosting bw file path, and link file path
                   eg. {'bw': {'Liver': 'path/to/Liver.bw', 'Kidney': 'path/to/Kidney.bw'},
                        'sashimiLinks': {'Liver': 'path/to/Liver.sashimi', 'Kidney': 'path/to/Kidney.sashimi'},
                        'unLinks': 'path/to/unprod.links',
                        'scale_link_heights': {'Liver': '0.5', 'Kidney': '0.5'}}
                        }
    nBins        : int. Number of bins for x axis of sashimiBigWig tracks

    return       : None
    """

    print(f"Generate ini files: {outIniFile} ...\n")
    with open(templateFile) as t:
        t_content = t.read()
        tt = Template(t_content)
        renderedContent = tt.render(
            group1=list(trackFiles["bw"].keys())[0],
            group1_bw=list(trackFiles["bw"].values())[0].split("/")[
                -1
            ],  # only file name
            group1_lnk=list(trackFiles["sashimiLinks"].values())[0].split("/")[-1],
            group1_lnk_ht_scale=list(trackFiles["scale_link_heights"].values())[0],
            n_bins=nBins,
            group2=list(trackFiles["bw"].keys())[1],
            group2_bw=list(trackFiles["bw"].values())[1].split("/")[-1],
            group2_lnk=list(trackFiles["sashimiLinks"].values())[1].split("/")[-1],
            group2_lnk_ht_scale=list(trackFiles["scale_link_heights"].values())[1],
            up_link_file=trackFiles["upLinks"].split("/")[-1],
        )
    with open(outIniFile, "w") as f:
        f.write(renderedContent + "\n")


def writePlotShell(templateFile, outShellFile, iniFile, plotFile, plotRange, plotTitle):
    """Generate shell script to plot sashimi plots
    templateFile : str. Path to template shell script
    outShellFile : str. Path to write shell script
    iniFile      : str. Path to ini file
    plotFile     : str. Path to output plot file
    plotRange    : tuple. Genomic range, BED like 0-based coordinates,
                   eg. chr1:25101-27101
    plotTitle    : str. Title of the plot, eg: 'Liver v Kidney: clu_261_+'

    """

    print(f"Generate shell script: {outShellFile} ...\n")
    with open(templateFile) as t:
        t_content = t.read()
        tt = Template(t_content)
        renderedContent = tt.render(
            iniFile=iniFile, plotFile=plotFile, plotRange=plotRange, plotTitle=plotTitle
        )

    with open(outShellFile, "w") as f:
        f.write(renderedContent + "\n")


def interactive_run():

    # for debugging

    OutDirPrefix = "test-prepSashimi"

    Contrast = "Brain-Cortex_v_Muscle-Skeletal"
    Groups = Contrast.split("_v_")

    # available bigwig files for given groups
    BwDirs = [f"resources/GTEx/BigWig/{x}" for x in Groups]
    BwFiles = {x: getBigWigFiles(y) for x, y in zip(Groups, BwDirs)}

    # sampels used in diff. splicing
    DsSampleGroupFile = f"results/ds/GTEx/{Contrast}/ds_sample_group.txt"
    DsSamples = getSamplesByGroup(DsSampleGroupFile, Groups)

    # plotable bigwig files are intersect of samples in ds and available bigwig files
    BwPlot = {x: getPlottableBigWig(BwFiles[x], DsSamples[x]) for x in Groups}

    # leafcutter output introns - used to just get intron type,
    IntronCountsFile = (
        f"results/ds/GTEx/{Contrast}/ds_perind_numers.counts.noise_by_intron.gz"
    )
    Introns = readIntronTable(IntronCountsFile)

    # PSI values from ds effect sizes
    DsEffectFile = f"results/ds/GTEx/{Contrast}/ds_effect_sizes.txt"
    EffectSizes = pd.read_csv(DsEffectFile, sep="\t")

    # plot target - use cluster id to filter introns
    PlotIntron = "chr1:15651895:15656420:clu_261_+"

    PlotClu = PlotIntron.split(":")[-1]
    PlotRange = getPlotRange(PlotClu, Introns, minRange=20000)

    # get links for sashimi track and unprod. links
    links = getPSI(EffectSizes, PlotClu, Introns, Groups)

    TrackDict = {
        "bw": {},
        "sashimiLinks": {},
        "upLinks": "",
        "scale_link_heights": {},
    }

    # write sashimi links
    for grp in links["sashimiLinks"]:
        outSashimiLinks = f"{OutDirPrefix}/{Contrast}/{grp}_{PlotClu}_{':'.join([str(x) for x in PlotRange])}.sashimi.links"
        TrackDict["sashimiLinks"][grp] = outSashimiLinks
        try:
            links["sashimiLinks"][grp].to_csv(
                outSashimiLinks, sep="\t", header=False, index=False
            )
        except OSError:
            os.makedirs(os.path.dirname(outSashimiLinks))
            links["sashimiLinks"][grp].to_csv(
                outSashimiLinks, sep="\t", header=False, index=False
            )

    # write unproductive links
    outUnprodLinks = f"{OutDirPrefix}/{Contrast}/{PlotClu}_{':'.join([str(x) for x in PlotRange])}.unproductive.links"
    TrackDict["upLinks"] = outUnprodLinks
    try:
        links["upLinks"].to_csv(outUnprodLinks, sep="\t", header=False, index=False)
    except OSError:
        os.makedirs(os.path.dirname(outUnprodLinks))
        links["upLinks"].to_csv(outUnprodLinks, sep="\t", header=False, index=False)

    # average bigwig files
    for grp in BwPlot:
        bw_files = BwPlot[grp]
        avgBw = avgBigwig(bw_files.values(), PlotRange)
        TrackDict["scale_link_heights"][grp] = np.round(
            np.mean(avgBw["values"]) / np.max(avgBw["values"]), 5
        )
        outBw = f"{OutDirPrefix}/{Contrast}/{grp}_{PlotClu}_{':'.join([str(x) for x in PlotRange])}.bw"
        TrackDict["bw"][grp] = outBw
        writeBigwig(outBw, avgBw, PlotRange)

    # write ini file
    TemplateFile = "test-prepSashimi/template-dge-ds.ini"
    OutIniFile = f"{OutDirPrefix}/{Contrast}/{PlotClu}_{':'.join([str(x) for x in PlotRange])}.ini"
    writeIni(TemplateFile, OutIniFile, TrackDict)

    writePlotShell(
        templateFile="test-prepSashimi/plot-sashimi-template.sh",
        outShellFile=f"{OutDirPrefix}/{Contrast}/{PlotClu}_{':'.join([str(x) for x in PlotRange])}.sh",
        iniFile=OutIniFile.split("/")[-1],
        plotFile=f"{PlotClu}_{':'.join([str(x) for x in PlotRange])}.pdf",
        plotRange=f"{PlotRange[0]}:{str(PlotRange[1])}-{str(PlotRange[2])}",
        plotTitle=f'"{Contrast}: {PlotClu}"',
    )


def main(args):
    """
    args: argpase object
    """

    (
        contrast,
        outDir,
        bwPrefix,
        dsSampleGroupFile,
        dsEffectFile,
        intronCountsFile,
        plotIntron,
        iniTemplate,
        plotShellTemplate,
    ) = (
        args.contrast,
        args.outDir,
        args.bwPrefix,
        args.dsSampleGroupFile,
        args.dsEffectFile,
        args.intronCountsFile,
        args.plotIntron,
        args.iniTemplate,
        args.plotShellTemplate,
    )

    Groups = contrast.split("_v_")

    # available bigwig files for given groups
    BwDirs = [f"{bwPrefix}/{x}" for x in Groups]
    BwFiles = {x: getBigWigFiles(y) for x, y in zip(Groups, BwDirs)}

    # samples in diff. splicing
    DsSamples = getSamplesByGroup(dsSampleGroupFile, Groups)

    # plotable bigwig files are intersect of samples in ds and available bigwig files
    BwPlot = {x: getPlottableBigWig(BwFiles[x], DsSamples[x]) for x in Groups}

    # leafcutter output introns - used to just get intron type
    Introns = readIntronTable(intronCountsFile)

    # PSI values from ds effect sizes
    EffectSizes = pd.read_csv(dsEffectFile, sep="\t")

    # ----------------- Plot target -----------------#
    PlotClu = plotIntron.split(":")[
        -1
    ]  # get clu id from intron, eg. chr1:15651895:15656420:clu_261_+
    PlotRange = getPlotRange(PlotClu, Introns, minRange=20000)

    # get links for sashimi track and unprod. links
    links = getPSI(EffectSizes, PlotClu, Introns, Groups)

    # store bw files and sashimi links and unprod. links
    TrackDict = {
        "bw": {},
        "sashimiLinks": {},
        "upLinks": "",
        "scale_link_heights": {},
    }

    # write sashimi links
    for grp in links["sashimiLinks"]:
        outSashimiLinks = f"{outDir}/{grp}_{PlotClu}_{':'.join([str(x) for x in PlotRange])}.sashimi.links"
        TrackDict["sashimiLinks"][grp] = outSashimiLinks
        try:
            links["sashimiLinks"][grp].to_csv(
                outSashimiLinks, sep="\t", header=False, index=False
            )
        except OSError:
            os.makedirs(os.path.dirname(outSashimiLinks))
            links["sashimiLinks"][grp].to_csv(
                outSashimiLinks, sep="\t", header=False, index=False
            )

    # write unproductive links
    outUnprodLinks = (
        f"{outDir}/{PlotClu}_{':'.join([str(x) for x in PlotRange])}.unproductive.links"
    )
    TrackDict["upLinks"] = outUnprodLinks
    try:
        links["upLinks"].to_csv(outUnprodLinks, sep="\t", header=False, index=False)
    except OSError:
        os.makedirs(os.path.dirname(outUnprodLinks))
        links["upLinks"].to_csv(outUnprodLinks, sep="\t", header=False, index=False)

    # average bigwig files
    for grp in BwPlot:
        bw_files = BwPlot[grp]
        avgBw = avgBigwig(bw_files.values(), PlotRange)
        outBw = f"{outDir}/{grp}_{PlotClu}_{':'.join([str(x) for x in PlotRange])}.bw"
        TrackDict["bw"][grp] = outBw
        TrackDict["scale_link_heights"][grp] = np.round(
            np.mean(avgBw["values"]) / np.max(avgBw["values"]), 5
        )
        writeBigwig(outBw, avgBw, PlotRange)

    OutIniFile = f"{outDir}/{PlotClu}_{':'.join([str(x) for x in PlotRange])}.ini"
    writeIni(iniTemplate, OutIniFile, TrackDict)

    writePlotShell(
        templateFile=plotShellTemplate,
        outShellFile=f"{outDir}/{PlotClu}_{':'.join([str(x) for x in PlotRange])}.sh",
        iniFile=OutIniFile.split("/")[-1],
        plotFile=f"{PlotClu}_{':'.join([str(x) for x in PlotRange])}.pdf",
        plotRange=f"{PlotRange[0]}:{str(PlotRange[1])}-{str(PlotRange[2])}",
        plotTitle=f'"{contrast}\n{PlotClu}:{PlotRange[0]}:{str(PlotRange[1])}-{str(PlotRange[2])}"',
    )


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description="Prepare dataset for sashimi plots")
    parser.add_argument(
        "--contrast",
        type=str,
        help='contrast in differential splicing, eg. "Brain-Cortex_v_Muscle-Skeletal"',
    )
    parser.add_argument(
        "--outDir", type=str, help='output directory prefix, eg. "test-prepSashimi"'
    )
    parser.add_argument(
        "--bwPrefix",
        type=str,
        help='prefix of directory containing bigwig files, eg. "resources/GTEx/BigWig"',
    )
    parser.add_argument(
        "--dsSampleGroupFile",
        type=str,
        help='sample group file from differential splicing, eg. "results/ds/GTEx/Brain-Cortex_v_Muscle-Skeletal/ds_sample_group.txt"',
    )
    parser.add_argument(
        "--dsEffectFile",
        type=str,
        help='effect size file from differential splicing, eg. "results/ds/GTEx/Brain-Cortex_v_Muscle-Skeletal/ds_effect_sizes.txt"',
    )
    parser.add_argument(
        "--intronCountsFile",
        type=str,
        help='intron counts file, eg. "results/ds/GTEx/Brain-Cortex_v_Muscle-Skeletal/ds_perind_numers.counts.noise_by_intron.gz"',
    )
    parser.add_argument(
        "--plotIntron",
        type=str,
        help='intron to plot, eg. "chr1:15651895:15656420:clu_261_+"',
    )
    parser.add_argument(
        "--iniTemplate",
        type=str,
        help="template ini file for pyGenomeTracks",
    )
    parser.add_argument(
        "--plotShellTemplate",
        type=str,
        help="template shell script for plotting sashimi plots",
    )

    args = parser.parse_args()

    main(args)
