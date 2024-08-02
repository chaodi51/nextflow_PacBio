#!/usr/bin/env python

""" Isolate clipped flanks of transposon alignments and output as fastq """

import argparse

import pysam

parser = argparse.ArgumentParser()
parser.add_argument('--bam')
parser.add_argument('--fastq')
parser.add_argument('--tsv')
args = parser.parse_args()

pysam.index(args.bam)
bam = pysam.AlignmentFile(args.bam)
reflen = max(bam.lengths)
alignments = bam.fetch()

with open(args.tsv, 'w') as tsv:
    tsv.write('\t'.join(['read', 'transposon', 'tn_start', 'tn_end']) + '\n')

with open(args.fastq, 'w') as fq:
    pass

for align in alignments:

    if align.is_supplementary:
        continue

    readname = align.query_name
    refname = align.reference_name
    sequence = align.query_sequence
    qualities = quality = ''.join(map(lambda x: chr( x+33 ), align.query_qualities))
    refpos = align.get_reference_positions()
    align_start = min(refpos)
    align_end = max(refpos)
    
    with open(args.tsv, 'a') as tsv:
        tsv.write('\t'.join([readname, refname, str(align_start), str(align_end)]) + '\n')

    cigartuples = align.cigartuples
    for index, (op, run) in enumerate(cigartuples):
        if op == 4:
            if index == 0:
                fqname = readname + '/left'
                fqseq = sequence[:run]
                fqqual = qualities[:run]
            else:
                fqname = readname + '/right'
                fqseq = sequence[-run:]
                fqqual = qualities[-run:]
            with open(args.fastq, 'a') as fq:
                fq.write('\n'.join(['@'+fqname, fqseq, '+', fqqual]) + '\n')
