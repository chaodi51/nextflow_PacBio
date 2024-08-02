#!/usr/bin/env python3

import argparse
import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument('--alvis')
parser.add_argument('--workflow-id')
parser.add_argument('--output')

args = parser.parse_args()

chimera_table = pd.read_csv(args.alvis, sep='\t', names=["read_id", "mid_pos", "ref1", "ref2", "read_p1_st", "read_p1_end",
    "read_p2_st", "read_p2_end", "align_p1_st", "align_p1_end", "align_p2_st", "align_p2_end"])
chimera_table = chimera_table[~(chimera_table['ref1'] == chimera_table['ref2'])] 
chimera_table['workflow_id'] = args.workflow_id
chimera_table.set_index("read_id", inplace=True)
chimera_table.to_csv(args.output, sep='\t')
