#!/usr/bin/env python


import argparse
import os
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="File containing contigs|scaffolds to extract sequences from")
parser.add_argument("-b", "--binning", help="File containing binning")
parser.add_argument("-o", "--output", help="Directory to save resulting bins in")
args = parser.parse_args()

contigs = list(SeqIO.parse(args.input, "fasta"))
binning_df = pd.read_table(args.binning, header=None, names=["SEQUENCEID", "BINID"], comment="@")
binning_df = binning_df.groupby("BINID").aggregate(set)

for idx, row in binning_df.iterrows():
    records_to_save = filter(lambda rec: rec.id in row.SEQUENCEID, contigs)
    with open(os.path.join(args.output, f"{idx}.fasta"), "w") as file:
        SeqIO.write(records_to_save, file, "fasta")

