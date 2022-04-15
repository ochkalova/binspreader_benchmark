#!/usr/bin/env python


import argparse
import os
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument("--binning_table", default=None, help="TSV file containing 2 fields: SEQUENCEID and BINID")
parser.add_argument("--binning_files", default=None, nargs="*", help="List of binning fasta files")
args = parser.parse_args()


if args.binning_table:
    current_bin = 1
    old_new_names_map = {}
    binning_df = pd.read_table(args.binning_table, header=None, names=["SEQUENCEID", "BINID"], comment="@")
    for idx, row in binning_df.iterrows():
        if row.BINID not in old_new_names_map:
            old_new_names_map[row.BINID] = f"bin_{current_bin}"
            current_bin += 1
        binning_df.loc[idx, "BINID"] = old_new_names_map[row.BINID]

    binning_df.to_csv(args.binning_table, index=None, sep="\t", header=None)

if args.binning_files:
    current_bin = 1
    for bin_file in args.binning_files:
        os.rename(bin_file, os.path.join(os.path.dirname(bin_file), f"bin_{current_bin}.fasta"))
        current_bin += 1

