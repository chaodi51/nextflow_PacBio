#!/usr/bin/env python3

import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--bed')
parser.add_argument('--uniqueid')
args = parser.parse_args()

colnames = ['reference', 'start', 'end', 'firstfrac', 'totalstr', 'totalcount', 'methyl', 'nonmethyl', 'methylfrac'] 
df = pd.read_csv(args.bed, sep='\t', names=colnames, usecols=['reference', 'start', 'totalcount', 'methylfrac'])
df['workflow_id'] = args.uniqueid

df.to_csv('methylation.tsv', sep='\t', index=False)

