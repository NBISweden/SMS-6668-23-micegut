#!/usr/bin/env python

import pandas as pd
from glob import glob
from argparse import ArgumentParser
import sys


def main(args):
    samples = {}
    files = glob(args.indir + "/*/assembly/contig_stats/postfilter_coverage_stats.txt")
    for f in files:
        sample = f.split("/")[-4]
        df = pd.read_csv(f, sep="\t", index_col=0)
        mapped = df["Plus_reads"].sum() + df["Minus_reads"].sum()
        samples[sample] = mapped
    df = pd.DataFrame(samples, index=["Mapped_reads"])
    df.index.name = "Sample"
    df.T.to_csv(sys.stdout, sep="\t")

if __name__ == "__main__":
    parser = ArgumentParser(
        description="""
This script searches for files matching pattern '<indir>/*/assembly/contig_stats/postfilter_coverage_stats.txt' and outputs total number of reads mapped.
        """
    )
    parser.add_argument("indir", help="Path to directory containing samples")
    args = parser.parse_args()
    main(args)