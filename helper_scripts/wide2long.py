#!/usr/bin/env python

import argparse
import io
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument("binning_file")
parser.add_argument("--add-header", default=False, action="store_true", help="Add CAMI format header")
parser.add_argument("--sample-id", type=str, default="SAMPLEID", help="Sample ID for CAMI format header")
parser.add_argument("--skip_zero_bin", default=False, action="store_true", help="Ignore bin 0")
args = parser.parse_args()

output_file_buffer = io.StringIO()

if args.add_header:
    print(f"@Version:0.9.0\n@SampleID:{args.sample_id}\n@@SEQUENCEID\tBINID", file=output_file_buffer)


contig_bin_pairs = []
with open(args.binning_file) as file_in:
    for line in file_in:
        if line[0] in {"@", "#"}:   # Skip headers
            print(line.strip(), file=output_file_buffer)
            continue
        contig_name, *bins = line.strip().split()
        if bins[0] == '0' and args.skip_zero_bin:
            continue
        for bin_ in bins:
            contig_bin_pairs.append((contig_name, bin_))

data = pd.DataFrame.from_records(contig_bin_pairs)
data.to_csv(output_file_buffer, header=None, index=None, sep="\t")
output_file_buffer.seek(0)
print(output_file_buffer.read().strip())
