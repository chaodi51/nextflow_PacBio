#!/usr/bin/env python3

import argparse
import glob
import json
import os
import requests
import subprocess
import sys
import warnings

import pandas as pd

from time import sleep

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))

def replace_in_file(infile, outfile, replace_dict):
    """ In infile, replace instances of one or more strings and write
            the resulting modified file to outfile
        replace_dict matches strings to be replaced with their replacements
    """
    with open(infile, 'r') as inlines:
        with open(outfile, 'w') as outlines:
            for line in inlines:
                newline = line
                for instring, outstring in replace_dict.items():
                    newline = newline.replace(instring, outstring)
                outlines.write(newline)

def addquotes(string):
    """ Put double quotes around a string """
    return '"' + string + '"'

def get_deployment_url(stdout):
    """ Extract content URL from deployment script standard output
        Returns None if URL cannot be found
    """
    for line in stdout.split('\n'):
        if line.startswith('Document successfully deployed'):
            url = line.rstrip().split(' ')[-1]
            return url
    return None

def set_tags(server_tags, wanted_tags, api_url, headers):
    """ Set wanted_tags on published content 
        api_url is the API url, not that at which the content is viewed
        server_tags are all tags known to server
        headers contains Connect credentials
    """
    for tag in server_tags:
        if tag['name'] in wanted_tags:
            response = requests.post(api_url+'/tags', headers=headers, json={'tag_id': tag['id']})
    response = requests.get(api_url+'/tags', headers=headers).json()
    posted_tags = [t['name'] for t in response]
    for tag in wanted_tags:
        if tag not in posted_tags:
            warnings.warn('Intended tag %s not added' % tag)

def set_access_all(url, domain, api_endpoint, headers):
    """ Set access to 'all' for content published to a url
        Returns guid for that content
    """
    access_all = {'access_type': 'all'}
    guid = url.split('/')[-2]
    api_url = ''.join([domain, api_endpoint, 'experimental/content/', guid])
    response = requests.post(api_url, headers=headers, json=access_all)
    return guid

# params.outdir = "s3://sparktx-nextflow-pipeline/projects/nf_pacbio/outputs/${params.unique_id}/"
# params.outdir = "${projectDir}/debug"
# params.datalakedir = "s3://sparkds-datalake-groupdropin-bioinformatics/nf_pacbio/${params.unique_id}/"

DEST_S3_BUCKET = 's3://sparkds-datalake-groupdropin-bioinformatics/nf_pacbio/'
DOMAIN = 'https://connect.sparkds.io'
API_ENDPOINT = '/__api__/v1/'
API_USER = os.environ.get("API_USER")
API_KEY = os.environ.get("API_KEY")

HEADER = {'Authorization': 'Key ' + API_KEY}

parser = argparse.ArgumentParser()
parser.add_argument('--version', default=None, help='Workflow version of run to be published, e.g. 1.0.6')
parser.add_argument('--sessionid', help='Tower sessionId for the run to be published')
parser.add_argument('--report-type', default='full', help='Type of report to be published')
parser.add_argument('--report-template', default=SCRIPT_DIR+'/../templates/template_nextflow_report.Rmd', help='Rmd template for sample report')
parser.add_argument('--toc-template', default=SCRIPT_DIR+'/../templates/template_nextflow_aggregate.Rmd', help='Rmd template for table of contents')
parser.add_argument('--report-script', default=SCRIPT_DIR+'/../Rscripts/deploy_nextflow_sample.R', help='R script to deploy sample report')
parser.add_argument('--toc-script', default=SCRIPT_DIR+'/../Rscripts/deploy_nextflow_agg.R', help='R script to deploy table of contents')
parser.add_argument('--debug', default=False, action='store_true', help='Set this to turn on print statements for debugging')
args = parser.parse_args()

VERSION = args.version
TABLE_OF_CONTENTS = 'nextflow_deployment_toc_%s.tsv' % args.sessionid

deployment_info = {}

metadata_loc = '%sv%s/%s/metadata.tsv' % (DEST_S3_BUCKET, args.version, args.sessionid)
aws_shared_loc = '%sv%s/' % (DEST_S3_BUCKET, args.version)
local_metadata = args.sessionid+'_metadata.tsv'

metadata_proc = subprocess.run(['aws', 's3', 'cp', metadata_loc, local_metadata, '--sse', 'AES256', '--acl', 'bucket-owner-full-control'], capture_output=True)
if args.debug:
    print(metadata_proc.stdout.decode())
    print(metadata_proc.stderr.decode())
metadata = pd.read_csv(local_metadata, sep='\t')
sample = metadata['sample_id'][0]

fasta_loc = metadata['ref_path'][0] + '/' + os.path.basename(metadata['primary_ref'][0]) + '.fasta'
aws_fasta_dir = '%srefs/' % (DEST_S3_BUCKET)
fasta_proc = subprocess.run(['aws', 's3', 'cp', fasta_loc, aws_fasta_dir, '--sse', 'AES256', '--acl', 'bucket-owner-full-control'], capture_output=True)
if args.debug:
    print(fasta_proc.stdout.decode())
    print(fasta_proc.stderr.decode())

all_tags = requests.get(''.join([DOMAIN, API_ENDPOINT, 'tags']), headers=HEADER).json()
new_tags = ('TDO', 'NGS/Proteomics', metadata['program'][0], 'DNA-Sequencing', 'Adeno-associated Virus 2')

run = metadata['run'][0]
primary = metadata['primary_ref'][0]
unique_id = metadata['workflow_id'][0]
version = metadata['version'][0]

replacement_strings = {'VERSION_HERE': addquotes(VERSION), 
                       'SESSIONID_HERE': addquotes(args.sessionid)}

sample_rmd = os.path.join(os.path.split(os.path.abspath(args.report_template))[0], 'make_report.Rmd')
replace_in_file(args.report_template, sample_rmd, replacement_strings)
deployproc = subprocess.run(['Rscript', args.report_script, 
        '--run', str(unique_id), "--rmdfile", sample_rmd], capture_output=True)

if args.debug:
    print(deployproc.stdout.decode())
    print(deployproc.stderr.decode())
url = get_deployment_url(deployproc.stdout.decode())
if not url:
    raise RuntimeError('Deployment failed for %s %s' % (run, sample))

guid = set_access_all(url, DOMAIN, API_ENDPOINT, HEADER)
deployment_info['run'] = run
deployment_info['sample'] = sample
deployment_info['url'] = url
deployment_info['guid'] = guid
deployment_info['tower_id'] = unique_id
deployment_info['target'] = primary
deployment_info['project'] = metadata['project'][0]
deployment_info['group'] = metadata['group'][0]
deployment_info['department'] = metadata['department'][0]
deployment_info['requester'] = metadata['requester'][0]
deployment_info['date'] = metadata['date'][0]
deployment_df = pd.DataFrame([deployment_info])

api_url = ''.join([DOMAIN, API_ENDPOINT, 'content/', guid])
set_tags(all_tags, new_tags, api_url, HEADER)

deployment_df['report'] = deployment_df.apply(lambda row: "<a href=" + row["url"] + ">sample_report</a>", axis=1)
deployment_df.drop(['url', 'guid'], axis=1, inplace=True)
columnorder = ['report', 'tower_id', 'run', 'sample', 'target', 'project', 'group', 'department', 'requester', 'date']
deployment_df = deployment_df[columnorder]
deployment_df.to_csv(TABLE_OF_CONTENTS, sep='\t')
toc_proc = subprocess.run(['aws', 's3', 'cp', TABLE_OF_CONTENTS, aws_shared_loc+'rsconnect/', '--sse', 'AES256', '--acl', 'bucket-owner-full-control'], capture_output=True)
if args.debug:
    print(toc_proc.stdout.decode())
    print(toc_proc.stderr.decode())

toc_rmd = os.path.join(os.path.split(os.path.abspath(args.report_template))[0], 'make_aggregate.Rmd')
replace_in_file(args.toc_template, toc_rmd, replacement_strings)

while True:
    aggproc = subprocess.run(['Rscript', args.toc_script, '--version', version.replace('.', '_'), "--rmdfile", toc_rmd], capture_output=True)
    if 'superseded' in aggproc.stderr.decode() or 'superseded' in aggproc.stdout.decode():
        sleep(30)
    else:
        break
if args.debug:
    print(aggproc.stdout.decode())
    print(aggproc.stderr.decode())
aggurl = get_deployment_url(aggproc.stdout.decode())
if not aggurl:
    raise RuntimeError('Deployment failed for table of contents')
guid = set_access_all(aggurl, DOMAIN, API_ENDPOINT, HEADER)
api_url = ''.join([DOMAIN, API_ENDPOINT, 'content/', guid])
set_tags(all_tags, new_tags, api_url, HEADER)
os.remove(local_metadata)
os.remove(TABLE_OF_CONTENTS)
print('Table of contents deployed to %s' % aggurl)
