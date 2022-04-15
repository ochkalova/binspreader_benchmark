#!/usr/bin/env python3

import sys
import os
import gfapy
import argparse
from textwrap import wrap
from Bio.Seq import Seq

op = argparse.ArgumentParser(description=__doc__)
op.add_argument("filename")
op.add_argument("output")
opts = op.parse_args()


def inv_link(link):
    return (link[0], "+" if link[1] == "-" else "-")

def segment_name(gfa, name, orientation):
    sname = "EDGE_{}_length_{}".format(name,
                                       len(gfa.segment(name).sequence))
    DP = gfa.segment(name).try_get("DP")
    if DP:
        sname += "_cov_" + str(DP)
    return sname + ("'" if orientation == "-" else "")

gfa = gfapy.Gfa.from_file(opts.filename, vlevel=0)
try:
    links = {}

    for s in gfa.segments:
        links[(s.name, "+")] = []
        links[(s.name, '-')] = []
    
    for e in gfa.edges:
        source = (e.from_segment.name, e.from_orient)
        target = (e.to_segment.name, e.to_orient)
        links[source].append(target)
        links[inv_link(target)].append(inv_link(source))

    for (lfrom, lor), lto in links.items():
        inv = lor == "-"
        seq = Seq(gfa.segment(lfrom).sequence)
        hdr = ">{}".format(segment_name(gfa, lfrom, lor))
        if (len(lto) > 0):
            hdr += ":" + ",".join(map(lambda eo: segment_name(gfa, eo[0], eo[1]), lto))
        hdr += ";"
        print(hdr)

        print("\n".join(wrap(str(seq if not inv else seq.reverse_complement()), 60)))


except gfapy.Error as err:
    sys.stderr.write(str(err))
    sys.exit(1)
