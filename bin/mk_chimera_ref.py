#!/usr/bin/env python3
"""Rename chroms and combine fa files"""
import sys
from Bio import SeqIO


def add_fa(fout, prefix: str, fa: str):
    seq_iter = SeqIO.parse(fa, format="fasta")
    for r in seq_iter:
        name = ">" + prefix + "___" + r.name
        seq = str(r.seq)
        print(name, file=fout)
        print(seq, file=fout)


def main():
    host, aav, host_fasta, aav_fasta, out_fa = sys.argv[1:]
    ls = ((host, host_fasta), (aav, aav_fasta))
    with open(out_fa, "w") as fout:
        for prefix, fa in ls:
            add_fa(fout, prefix, fa)


main()
