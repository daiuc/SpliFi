#!/usr/bin/env python


__author__ = "Chao Dai"

import argparse
import gzip
import sys
from datetime import datetime

import pandas as pd

def pick_sjc_label(df):
    '''Pick one label for each intron based on priroty'''
    sjc_priority = {'PR': 1, 'UP': 2, 'NE':3} # PR > UP > NE
    df['priority'] = df['SJClass'].map(sjc_priority)
    chosen = df.sort_values(by = 'priority', ascending = True).iloc[0]
    return chosen.drop('priority')

def merge_discordant_logics(sjc_file: str):
    '''some junctions have multiple classifications. Use conservative approach
    to merge them.
    '''

    classifier = {
        # each bit represents:
        # is_GTFAnnotatedCoding?  is_GTFAnnotated?  is_LF2AnnotatedCoding?  is_ClosetoUTR?
        '0000': 'UP', # UnProductive,
        '0001': 'NE', # NEither productive nor unproductive
        '0010': 'PR', # PRoductive
        '0011': 'PR', # PRoductive
        '0100': 'UP', # UnProductive
        '0101': 'NE', # Neither Productive nor UnProductive
        '0110': 'PR', # PRoductive
        '0111': 'PR', # PRodutive
        '1000': 'PR', # PRoductive
        '1001': 'PR', # PRoductive
        '1010': 'PR', # PRoductive
        '1011': 'PR', # PRoductive
        '1100': 'PR', # PRoductive
        '1101': 'PR', # PRoductive
        '1110': 'PR', # PRoductive
        '1111': 'PR'  # PRoductive
        }

    classifer_3bits = {
        # each bit represents:
        # is_GTFAnnotated?  is_LF2AnnotatedCoding?  is_ClosetoUTR?
        '000': 'UP', # UnProductive,
        '001': 'NE', # NEither productive nor unproductive
        '010': 'PR', # PRoductive
        '011': 'PR', # PRoductive
        '100': 'UP', # UnProductive
        '101': 'NE', # Neither Productive nor UnProductive
        '110': 'PR', # PRoductive
        '111': 'PR', # PRodutive
    }

    sjc = pd.read_csv(sjc_file, sep = "\t")

    # group dt; NOTE:ForwardSpliceJunctionClassifier has an extra Strand column, backward doesn't
    if 'Strand' in sjc.columns:
        sjc = sjc[['Gene_name', 'Intron_coord', 'Strand', 'GencodePC', 'Annot', 'Coding', 'UTR']]

        # convert Annotation, Coding, UTR status to SJ categories
        sjc['SJClass'] = sjc[['GencodePC', 'Annot', 'Coding', 'UTR'
                              ]].astype(int).astype(str).agg(''.join, axis=1).map(classifier)

        # if multiple classifications, take the one with highest priority
        sjc = sjc.groupby(['Intron_coord', 'Strand']).apply(pick_sjc_label).reset_index(drop=True)
        sjc = sjc[['Intron_coord', 'Strand', 'SJClass', 'Gene_name']]

        # convert df to dict
        sjc = sjc.set_index(['Intron_coord', 'Strand']).to_dict(orient='index')
    else:
        sjc = sjc[['Gene_name', 'Intron_coord', 'Annot', 'Coding', 'UTR']]
        
        # convert Annotation, Coding, UTR status to SJ categories
        sjc['SJClass'] = sjc[['Annot', 'Coding', 'UTR'
                              ]].astype(int).astype(str).agg(''.join, axis=1).map(classifer_3bits)
        
        # if multiple classifications, take the one with highest priority
        sjc = sjc.groupby(['Intron_coord', 'Strand']).apply(pick_sjc_label).reset_index(drop=True)
        sjc = sjc[['Intron_coord', 'Strand', 'SJClass', 'Gene_name']]

        # convert df to dict
        sjc = sjc.set_index('Intron_coord').to_dict(orient='index')

    sjc = {flatten_tuple(k): v for k, v in sjc.items()}
    sys.stderr.write(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Loaded junction annotation lookup.\n")

    # sjc is a dcitionary with:
    # - keys: intron coordinates, e.g. ('chr1', 1000, 2000, '+') or ('chr1', 1000, 2000) for backward
    # - values: a dictionary e.g. {'SJClass': 'UP', 'Gene_name': 'DNMBP'}
    return sjc

 
def boolean_to_bit(bool_vec):
    # Convert boolean vector to string of "1"s and "0"s
    bin_str = "".join(["1" if b else "0" for b in bool_vec])

    return bin_str


def flatten_tuple(t):
    # t: tuple like ('chr1:100-200', '+')' or str 'chr1:100-200'
    if isinstance(t, tuple):
        c, ab = t[0].split(":")
    elif isinstance(t, str):
        c, ab = t.split(":")
    a, b = ab.split("-")
    a, b = int(a), int(b)
    if isinstance(t, tuple):
        s = t[1]
        return((c, a, b, s)) # e.g. ('chr1', 100, 200, '+')
    else:
        return((c, a, b))


def annotate_noisy(inFile, outPrefix, sjcFile):
    """Annotate introns using previously processed junction classifications.

    Produces 3 files. Details in side-effects.

    Parameters:
    -----------
        inFile: str. e.g. wConst_perind.counts.gz
        outPrefix: str. output prefix (incl. path)
        sjcFile: str. inFile paring junction annotation file.

    Returns:
    --------
        return : no return. Use side-effects to write output files.

    Side-effects:
    -------------
        noisy numerators : output result file.
            Same as noisy counts, except here only reports the numerators.
            Also, used `*` to mark both start and end pseudo coordinates
            of noisy intron.

        noisy counts by intron : output file for diagnostics.
            Same count table as the output of sorted juncs count table,
            except that each intron is annotated with `F`, `N`, or `PN` to
            denote `putative_functional`, `noisy`, or `putative_noisy`.

    """

    sys.stderr.write(
        f"\n##{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Annotating introns with custom-classified annotations {sjcFile}...\n"
    )

    sjc = merge_discordant_logics(sjcFile)

    noisydiag = f"{outPrefix}.constcounts.annotated.gz"
    numersdiag = f"{outPrefix}.numers.annotated.gz"
    foutdiag = gzip.open(noisydiag, "wt")
    foutdiagnumers = gzip.open(numersdiag, "wt")

    F = gzip.open(inFile)
    ln = F.readline().decode()
    foutdiag.write(ln)
    foutdiagnumers.write(ln)

    sys.stderr.write(f"##{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Annotating introns...\n")
    nline = 1
    for ln in F:
        if nline % 5000 == 0:
            sys.stderr.write(f"##{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} {nline} lines annotated.\n")
        if type(ln) == bytes:
            ln = ln.decode("utf-8")  # convert bytes to string
        ln = ln.split()
        intron = ln[0]
        chrom, s, e, clu = intron.split(":")  # chr, start, end, clu_1_+
        strand = clu.split("_")[-1]
        intronid = chrom, int(s), int(e), strand

        reads = [x.split("/")[0] for x in ln[1:]]  # numerators

        # annotate using custom classification
        if intronid in sjc:
            classification = sjc[intronid]["SJClass"]
        else:
            classification = "IN"  # IN: INtergenic

        # add class flag and write to *_perind.noise_by_intron.gz, eg: 'chr1:825552:829002:clu_1_+:F 1/14 0/25 1/33 1/14 1/33'
        foutdiag.write(intron + f":{classification}" + " " + " ".join(ln[1:]) + "\n")
        foutdiagnumers.write(
            intron
            + f":{classification}"
            + " "
            + " ".join([x for x in reads])
            + "\n"
        )

        nline += 1

    sys.stderr.write(f"##{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} {nline} lines annotated.\n")
    sys.stderr.write(f"##{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Annotation done.\n")
    foutdiag.close()


def main(options):

    inFile, outPrefix, sjcFile = options.inputFile, options.outPrefix, options.annotation
    annotate_noisy(inFile, outPrefix, sjcFile)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--inputFile",
        dest="inputFile",
        type=str,
        required=True,
        help="input perind.counts.gz file",
    )

    parser.add_argument(
        "-o",
        "--outPrefix",
        dest="outPrefix",
        type=str,
        required=True,
        help="Output prefix",
    )

    parser.add_argument(
        "-a",
        "--annotatio",
        dest="annotation",
        type=str,
        required=True,
        help="junction annotation file corresponding to the input file",
    )

    options = parser.parse_args()
    main(options)
