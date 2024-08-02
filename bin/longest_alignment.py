#!/usr/bin/env python3

import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--outfile')
parser.add_argument('--mergedfile')
parser.add_argument('--primaryref')
parser.add_argument('--countfile', nargs='*')
args = parser.parse_args()

all_counts = []
all_ref_names = []
for countfile in args.countfile:
    ref_name = '.'.join(countfile.split('/')[-1].split('.')[:-2])
    all_ref_names.append(ref_name)
    count_df = pd.read_csv(countfile, sep='\t', index_col='read')
    count_df.rename({'count': ref_name}, axis=1, inplace=True)
    all_counts.append(count_df)
merged_df = all_counts[0]
for df in all_counts[1:]:
    merged_df = merged_df.join(df, how='outer')

merged_df.fillna(0.0, inplace=True)

newcols = [args.primaryref] + [i for i in merged_df if i not in [args.primaryref]]
merged_df = merged_df[newcols]

merged_df['Longest'] = merged_df.idxmax(axis=1)
merged_df.to_csv(args.mergedfile, sep='\t')

counts_df = merged_df.Longest.value_counts()
for name in all_ref_names:
    if name not in counts_df.index:
        counts_df[name] = 0
counts_df.to_csv(args.outfile, sep='\t', header=False)
