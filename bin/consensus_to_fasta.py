#!/usr/bin/env python3

import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--tsv', help='consensus tsv')
parser.add_argument('--fasta', help='output fasta')
args = parser.parse_args()

df = pd.read_csv(args.tsv, sep='\t')

colnames = df.columns

max_consensus = 0
target_column = 'consensus'

"""
for name in colnames:
    if 'Consensus' in name:
        level = float(name.split(' ')[-1])
        if level > max_consensus:
            max_consensus = level
            target_column = name
if not target_column:
    raise IOError('No consensus column in input')
"""

consensus_vals = ['' if i == 'del' else i for i in list(df[target_column])]
consensus_seq = ''.join(consensus_vals)

with open(args.fasta, 'w') as fasta:
    fasta.write('>Consensus\n')
    fasta.write(consensus_seq)

