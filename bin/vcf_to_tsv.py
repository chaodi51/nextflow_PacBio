#!/usr/bin/env python

""" Transform VCF into a tabular format """

import argparse
import os

import pandas as pd

def parse_info(info_string):
    """ Break values out of INFO field and return
            as list, ints where possible or strings otherwise
    """
    parts = info_string.split(';')
    values = {}
    for field in parts:
        if '=' in field:
            key, value = field.split('=')
        else:
            key, value = field, 'True'
        values[key] = int_if_possible(value)
    return values

def int_if_possible(string):
    """ Return string cast to int if possible, otherwise
            original string
        Raises ValueError if input not of str type
    """
    if not isinstance(string, str):
        raise ValueError('Input not of str type')
    try:
        return int(string)
    except ValueError:
        return string

def verify_format_fields(record_format_fields_string, full_header):
    """ Verify that the FORMAT fields in a record (colon-
            delimited string) are the same as what appears in the
            header (last N fields, each preceded by 'FORMAT_')
        Returns bool
    """
    record_fields = record_format_fields_string.split(':')
    header_fields = full_header[-len(record_fields):]
    strip_header_fields = [f.split('_')[1] for f in header_fields]
    return record_fields == strip_header_fields 

parser = argparse.ArgumentParser()
#parser.add_argument('-s', '--sample', help='sample name')
#parser.add_argument('-r', '--run', help='run name')
parser.add_argument('-v', '--vcf', help='path to VCF')
#parser.add_argument('-t', '--vartype', help='SNP or INDEL')
parser.add_argument('-m', '--min-freq', type=float, help='minimum variant frequency to report call')
parser.add_argument('-f', '--frequency-called', help='name of expected allele frequence column')
parser.add_argument('--workflow_id')
parser.add_argument('-o', '--output', help='path to output tsv')

args = parser.parse_args()
if not os.path.isfile(args.vcf):
    raise IOError('Input VCF not found at %s' % args.vcf)
#if args.vartype not in ('SNP', 'INDEL'):
#    raise ValueError('vartype must be SNP or INDEL, not %s' % args.vartype)

#run = args.run
#sample = args.sample
vcf = args.vcf
#vartype = args.vartype
output = args.output

new_fields = []
format_exists = True
with open(vcf, 'r') as v:
    for line in v:
        skip = False
        # Get headers for INFO and FORMAT fields
        if line.startswith('##'):
            if not 'INFO' in line and not 'FORMAT' in line:
                continue
            field_name_info = line.lstrip('#').split(',')[0]
            fieldtype, _, name = field_name_info.split('=')
            new_fields.append('%s_%s' % (fieldtype, name))
        # Get standard headers
        elif line.startswith('#'):
            std_header = line.lstrip('#').rstrip().split('\t')
            info_pos = std_header.index('INFO')
            if not "FORMAT" in std_header:
                format_exists = False
            retained_fields = std_header[:info_pos]
            complete_header = retained_fields + new_fields
            df_dict = {field: [] for field in complete_header} 
        # Parse variant lines
        else:
            parts = line.rstrip().split('\t')
            retained_values = parts[:info_pos]
            info_values = parse_info(parts[info_pos])
            complete_info = [info_values[k[5:]] if k[5:] in info_values else '' for k in [i for i in complete_header if i.startswith('INFO')]]
            if format_exists:
                if not verify_format_fields(parts[info_pos+1], complete_header):
                    raise IOError('FORMAT inconsistency in %s' % vcf)
                format_values = parts[info_pos+2].split(':')
            else:
                format_values = []
            all_values = retained_values + complete_info + format_values
            if len(complete_header) != len(all_values):
                raise IOError('Number of fields does not match header')
            for key, value in zip(complete_header, all_values):
                if key == args.frequency_called and float(value) < args.min_freq:
                    skip = True
                    break
            if skip:
                continue
            for key, value in zip(complete_header, all_values):
                df_dict[key].append(value)

variant_df = pd.DataFrame(df_dict)
#variant_df.insert(0, 'type', vartype)
#variant_df.insert(0, 'sample', sample)
#if run:
#    variant_df.insert(0, 'run', run)
#variant_df = variant_df.astype({'FORMAT_RBQ': int, 'FORMAT_ABQ': int, 'FORMAT_RDF': int, 'FORMAT_ADF': int, 'FORMAT_RDR': int, 'FORMAT_ADR': int})
#variant_df['DELTA_BQ'] = variant_df['FORMAT_RBQ'] - variant_df['FORMAT_ABQ']
#variant_df['REF_FF'] = variant_df['FORMAT_RDF'] / (variant_df['FORMAT_RDF'] + variant_df['FORMAT_RDR'])
#variant_df['ALT_FF'] = variant_df['FORMAT_ADF'] / (variant_df['FORMAT_ADF'] + variant_df['FORMAT_ADR'])
if args.workflow_id:
    variant_df['workflow_id'] = args.workflow_id
variant_df.to_csv(output, sep='\t', index=False)
