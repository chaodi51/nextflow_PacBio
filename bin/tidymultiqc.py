#!/usr/bin/env python

""" A script for extracting QC metrics from multiqc json output and 
    formatting into a table
"""
import argparse
import json

import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--multiqc-json')
parser.add_argument('--run')
parser.add_argument('--workflow_id')
parser.add_argument('--output-tsv')

args = parser.parse_args()
multiqc_json = args.multiqc_json
run = args.run
output = args.output_tsv

with open(multiqc_json, 'r') as json_in:
    json_data = json.load(json_in)

general_data = json_data['report_general_stats_data']
fastqc_data = general_data[0]

fastqc_df = pd.DataFrame(fastqc_data).T

allmetrics = fastqc_df
allmetrics.reset_index(inplace=True)
allmetrics.rename({'index': 'sample_id'}, axis=1, inplace=True)
if run:
    allmetrics['run'] = run

allmetrics.dropna(how='any', inplace=True)

#if allmetrics['filtering_result_passed_filter_reads'].equals(allmetrics['after_filtering_total_reads']):
#    allmetrics.drop(columns=['after_filtering_total_reads'], inplace=True)
#else:
#    raise Exception('Read count metrics not redundant: %s' % tsv)

#allmetrics.index = [i[:-3] for i in allmetrics['sample_id']]
#allmetrics.drop(['sample_id', 'filtering_result_too_many_N_reads', 'filtering_result_too_short_reads', 'filtering_result_too_long_reads'], axis=1, inplace=True)
#allmetrics.index.rename('sample', inplace=True)

if args.workflow_id:
    allmetrics['workflow_id'] = args.workflow_id

allmetrics.reset_index(drop=True, inplace=True)
allmetrics.set_index("sample_id", inplace=True)

allmetrics.to_csv(output, sep='\t')
