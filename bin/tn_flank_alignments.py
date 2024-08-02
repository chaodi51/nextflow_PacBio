#!/usr/bin/env python

""" Isolate clipped flanks of transposon alignments and output as fastq """

import argparse

import pysam

parser = argparse.ArgumentParser()
parser.add_argument('--bam')
parser.add_argument('--output')
args = parser.parse_args()

pysam.index(args.bam)
bam = pysam.AlignmentFile(args.bam)
reflen = max(bam.lengths)
alignments = bam.fetch()

outputs = {'left': args.output + '_left.tsv', 'right': args.output + '_right.tsv'}
with open(outputs['left'], 'w') as left_out:
    left_out.write('\t'.join(['read', 'left_ref', 'left_start', 'left_end']) + '\n')
with open(outputs['right'], 'w') as right_out:
    right_out.write('\t'.join(['read', 'right_ref', 'right_start', 'right_end']) + '\n')

for align in alignments:
    if align.is_supplementary:
        continue
    readname = align.query_name
    chirality = readname.split('/')[-1]
    orig_readname = '/'.join(readname.split('/')[:-1])
    refname = align.reference_name
    refpos = align.get_reference_positions()
    align_start = min(refpos)
    align_end = max(refpos)
    with open(outputs[chirality], 'a') as output: 
        output.write('\t'.join([orig_readname, refname, str(align_start), str(align_end)]) + '\n')
