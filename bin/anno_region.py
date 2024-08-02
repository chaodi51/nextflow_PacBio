#!/usr/bin/env python3
'''
Add annotation to genome coverge file
'''
import re
import sys

bed, cov, out = sys.argv[1:]

pos_itr = {}
itr_border = []
with open(bed) as f:
    for line in f:
        if re.search("ITR", line):
            chr, st, end, anno, _, strand = line.strip().split("\t")
            itr_border = itr_border + [st, end]
            for i in range(int(st), int(end) + 1):
                pos_itr[i] = "ITR"

with open(cov) as fcov, open(out, "w") as fout:
    for line in fcov:
        pos = int(line.split("\t")[0])
        if pos in pos_itr:
            fout.write(line.strip() + "\t" + pos_itr[pos] + "\n")
        else:
            if pos > int(itr_border[1]) and pos < int(itr_border[2]):
                fout.write(line.strip() + "\tOthers" + "\n")