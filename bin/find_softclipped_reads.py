#!/usr/bin/env python

""" Identify reads with soft clipping > 10% of the reference length """

import argparse

import pysam

parser = argparse.ArgumentParser()
parser.add_argument('--bam')
args = parser.parse_args()

bam = pysam.AlignmentFile(args.bam)
reflen = max(bam.lengths)
alignments = bam.fetch()

for align in alignments:
    sc_count = 0
    cigartuples = align.cigartuples
    for ctuple in cigartuples:
        op, run = ctuple
        if op == 4:
            sc_count += run
    if sc_count > 0.1 * reflen:
        print(align.to_string())
