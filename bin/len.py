#!/usr/bin/env python3
import sys
import pysam

bamfile, len_tab = sys.argv[1:]

bam = pysam.AlignmentFile(bamfile)

with open(len_tab, 'w') as out:
    for b in bam:
        out.write(str(b.reference_length) + '\n')
