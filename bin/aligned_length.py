#!/usr/bin/env python3 

import argparse
import re
import sys

#parser = argparse.ArgumentParser()
#parser.add_argument('--reference')
#args = parser.parse_args()

print('\t'.join(['read','count']))
for line in sys.stdin:
    parts = line.split('\t')
    read = parts[0]
    cigar = parts[5]
    readlen = len(parts[9])
    Mstr = re.findall('(\d+)[M=]', cigar)
    aligned_count = sum([int(i) for i in Mstr])
    print('\t'.join([read, str(aligned_count)]))
