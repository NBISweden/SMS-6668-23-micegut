#!/usr/bin/env python

from argparse import ArgumentParser
import pandas as pd
from Bio.SeqIO import parse
import sys


def split_ranks(df, ranks):
    uniq_lin = df.lineage.unique()
    d = {}
    for lin in uniq_lin:
        d[lin] = dict(zip(ranks,lin.split(";")))
    return pd.DataFrame(d).T


def read_genecatalog(f):
    genecatalog = {}
    samples = {}
    for seq in parse(f, "fasta"):
        genecat_id = seq.id
        sample_part = seq.description.split(" ")[1]
        sample = sample_part.split("_")[0]
        contig = "_".join(sample_part.split("_")[0:2])
        try:
            samples[sample].append(contig)
        except KeyError:
            samples[sample] = [contig]
        genecatalog[genecat_id] = contig
    return genecatalog, samples


def main(args):
    genecatalog_faa = args.genecatalog_faa
    atlas_dir = args.atlas_dir

    # Read in gene catalog protein fasta file
    genecatalog = {}
    samples = {}
    genecatalog, samples = read_genecatalog(genecatalog_faa)
    catalogdf = pd.DataFrame.from_dict(genecatalog, orient="index", columns=["contig"])
    
    taxdf = pd.DataFrame()
    for sample in samples.keys():
        taxfile = f"{atlas_dir}/{sample}/annotation/mmseqs2/UniRef100.taxonomy.tsv"
        _taxdf = pd.read_csv(taxfile, sep="\t", header=None, index_col=0, usecols=[0,8], names=["contig","lineage"])
        _taxdf.fillna("uc", inplace=True)
        lineage2ranks = split_ranks(_taxdf, ["superkingdom","phylum","class","order","family","genus","species"])
        _taxdf = pd.merge(_taxdf, lineage2ranks, left_on="lineage", right_index=True)
        taxdf = pd.concat([taxdf, _taxdf])
        
    catalogdf = pd.merge(catalogdf, taxdf, left_on="contig", right_index=True, how="left")
    catalogdf.fillna("unknown", inplace=True)
    catalogdf.index.name="gene_id"
    catalogdf.loc[:, ["superkingdom","phylum","class","order","family","genus","species"]].to_csv(sys.stdout, sep="\t")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-g", "--genecatalog_faa", help="Gene catalog protein fasta file")
    parser.add_argument("-a", "--atlas_dir", help="Atlas output directory")
    args = parser.parse_args()
    main(args)