#!/usr/bin/env python3

import argparse
import pysam

parser = argparse.ArgumentParser()
parser.add_argument('--readnames')
parser.add_argument('--inputbam')
parser.add_argument('--outputreads')
args = parser.parse_args()

readnames = {}
with open(args.readnames, 'r') as rn:
    for line in rn:
        readnames[line.rstrip()] = True

pysam.index(args.inputbam)
bam = pysam.AlignmentFile(args.inputbam)
alignments = bam.fetch()

for align in alignments:
    if align.query_name in readnames:
        with open(args.outputreads, 'a') as output:
            output.write(align.to_string() + '\n')
