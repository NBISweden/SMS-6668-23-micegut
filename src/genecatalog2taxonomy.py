#!/usr/bin/env python

from argparse import ArgumentParser
import pandas as pd
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

def most_resolved(df):
    for rank in ["species","genus","family","order","class","phylum","superkingdom"]:
        cl = df.loc[(~df[rank].str.startswith("uc"))&(df[rank]!="unknown")]
        if cl.shape[0] > 0:
            assignment = cl.groupby(rank).size().sort_values().head(1).index[0]
            return df.loc[df[rank]==assignment].head(1)
    return df.head(1)
            


def main(args):
    orf_info = args.orf_info
    atlas_dir = args.atlas_dir

    # Read in gene catalog protein fasta file
    genecatalog = {}
    samples = {}
    sys.stderr.write("Reading gene catalog clustering info\n")
    genecatalog = pd.read_parquet(orf_info)
    sys.stderr.write("Renaming gene nrs\n")
    Ngenes = genecatalog.GeneNr.max()
    n_leading_zeros = len(str(Ngenes))
    number_format = f"Gene{{:0{n_leading_zeros}d}}"
    genecatalog.GeneNr = genecatalog.GeneNr.apply(number_format.format)
    sys.stderr.write("Creating contig column\n")
    genecatalog['contig'] = genecatalog[["Sample","ContigNr"]].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
    genecatalog.set_index("GeneNr", inplace=True)
    catalogdf = pd.DataFrame(genecatalog.loc[:, "contig"])
    
    taxdf = pd.DataFrame()
    sys.stderr.write(f"Reading taxonomic assignments for {len(genecatalog.Sample.unique())} samples\n")
    for sample in genecatalog.Sample.unique():
        taxfile = f"{atlas_dir}/{sample}/annotation/mmseqs2/UniRef100.taxonomy.tsv"
        _taxdf = pd.read_csv(taxfile, sep="\t", header=None, index_col=0, usecols=[0,8], names=["contig","lineage"])
        _taxdf.fillna("uc", inplace=True)
        lineage2ranks = split_ranks(_taxdf, ["superkingdom","phylum","class","order","family","genus","species"])
        _taxdf = pd.merge(_taxdf, lineage2ranks, left_on="lineage", right_index=True)
        taxdf = pd.concat([taxdf, _taxdf])
    
    sys.stderr.write("Merging gene catalog and taxonomic assignments\n")
    catalogdf = pd.merge(catalogdf, taxdf, left_on="contig", right_index=True, how="left")
    catalogdf.fillna("unknown", inplace=True)
    catalogdf.index.name="gene_id"
    sys.stderr.write("Selecting most resolved taxonomic assignment\n")
    catalogdf = catalogdf.groupby(level=0).apply(most_resolved)
    catalogdf.loc[:, ["superkingdom","phylum","class","order","family","genus","species"]].to_csv(sys.stdout, sep="\t")


if __name__ == "__main__":
    parser = ArgumentParser("""
python src/genecatalog2taxonomy.py --orf_info atlas/Genecatalog/gene_catalog.faa -a atlas > atlas/G
enecatalog/gene_catalog.UniRef100.taxonomy.tsv
                            """)
    parser.add_argument("--orf_info", help="Orf info file from clustering (e.g. atlas/Genecatalog/clustering/orf_info.parquet)")
    parser.add_argument("-a", "--atlas_dir", help="Atlas output directory")
    args = parser.parse_args()
    main(args)