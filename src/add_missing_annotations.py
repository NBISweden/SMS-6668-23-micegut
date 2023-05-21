#!/usr/bin/env python

from argparse import ArgumentParser
import os
from glob import glob


def get_mags(mag_dir):
    for mag os.listdir(mag_dir):
        if not os.path.exists(f"{mag_dir}/{mag}/rrnas.tsv"):
            with open(f"{mag_dir}/{mag}/rrnas.tsv", 'w') as fhout:
                fhout.write("\t".join(["scaffold","fasta","begin","end","strand","type","e-value","note"]))
        if not os.path.exists(f"{mag_dir}/{mag}/trnas.tsv"):
            with open(f"{mag_dir}/{mag}/trnas.tsv", "w") as fhout:
                fhout.write("\t".join(["fasta","Name","tRNA #","Begin","End","Type","Codon","Score","Note"]))
            

def main():
    parser = ArgumentParser()
    parser.add_argument("mag_dir", type=str,
                       help="Path to directory with MAG annotations (in subdirs)")
    args = parser.parse_args()
    add_missing(args.mag_dir)

    
if __name__ == "__main__":
    main()