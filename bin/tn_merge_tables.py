#!/usr/bin/env python3

import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--cov')
parser.add_argument('--uniqueid')
parser.add_argument('--flanks', nargs='*')
args = parser.parse_args()

cov_table = pd.read_csv(args.cov, sep='\t')
flank_table1 = pd.read_csv(args.flanks[0], sep='\t')
flank_table2 = pd.read_csv(args.flanks[1], sep='\t')

merged = cov_table.merge(flank_table1, how='outer').merge(flank_table2, how='outer').set_index('read')
merged['workflow_id'] = args.uniqueid

merged.to_csv('tn_metrics.tsv', sep='\t', na_rep='N/A')
