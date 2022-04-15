#!/usr/bin/env python

import argparse
from collections import defaultdict, Counter

def read_binning(binning_file):
    with open(binning_file, 'r') as read_handler:
        for line in read_handler:
            if line.startswith('@') or line.startswith('#'):
                continue
            line = line.strip().split('\t')
            binnacle_scaff = line[0]
            bin_id = line[1]
            yield binnacle_scaff, bin_id
            
            
def read_coordinates(coordinates_file):
    with open(coordinates_file, 'r') as read_handler:
        binnacle2orig = defaultdict(list)
        orig2binnacle = Counter()
        for line in read_handler:
            if line.startswith('@') or line.startswith('#'):
                continue
            line = line.strip().split('\t')
            binnacle_scaff = "Binnacle_Scaffold_%s" % line[0]
            original_scaff = line[2]
            binnacle2orig[binnacle_scaff].append(original_scaff)
            orig2binnacle[original_scaff] += 1 
    if orig2binnacle.most_common(1)[0][1] != 1:
        print("Error! Correspondence of scaffolds is ambiguous!")
        raise BaseException
    return binnacle2orig


def convert(coordinates, binning, output_file):
    with open(output_file, 'w') as write_handler:
        binnacle2orig = read_coordinates(coordinates)
       # write_handler.write("#CAMI Format for Binning\n@Version:0.9.0\n@SampleID:_SAMPLEID_\n@@SEQUENCEID\tBINID\n")
        for binnacle_scaff, bin_id in read_binning(binning):
            for scaffold in binnacle2orig[binnacle_scaff]:
                write_handler.write("%s\t%s\n" % (scaffold, bin_id))


def main():
    parser = argparse.ArgumentParser(description="Convert Binnacle scaffolds to original scaffolds in TSV file")
    parser.add_argument("--coordinates", help="Coords_After_Delinking.txt file from Binnacle output")
    parser.add_argument("--binning", help="Binning in TSV format from Binnacle")
    parser.add_argument("-o", "--output_file", required=False, help="Output tsv file")
    
    args = parser.parse_args()
    convert(args.coordinates, args.binning, args.output_file)


if __name__ == "__main__":
    main()
