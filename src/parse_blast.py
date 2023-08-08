#!/usr/bin/env python

import pandas as pd
from argparse import ArgumentParser
import sys


def sort_bitscore(df):
    return df.sort_values("bitscore", ascending=False)


def main(args):
    header = ["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"]
    df = pd.read_csv(args.blastfile, sep="\t", header=None, names=header)
    df_filt = df.loc[(df.evalue<args.evalue)&(df.pident>=args.ident)]
    df_filt_sort = df_filt.groupby(level=0).apply(sort_bitscore)
    df_filt_best = df_filt_sort.groupby("qseqid").head(1).set_index("qseqid")
    df_filt_best.to_csv(sys.stdout, sep="\t")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("blastfile", type=str,
                        help="Blast output")
    parser.add_argument("-e", "--evalue", type = float, default = 10,
                        help="Evalue cutoff")
    parser.add_argument("-i", "--ident", type = int, default = 90,
                        help="Percent identity cutoff")
    args = parser.parse_args()
    main(args)
