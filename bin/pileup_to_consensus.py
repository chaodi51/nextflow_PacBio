#!/usr/bin/env python

import argparse
import pandas as pd
import re
import sys

class Variant(object):
    def __init__(self, pos, ref, alt, refdepth, altdepth):
        self.pos = int(pos)
        self.ref = ref
        self.alt = alt
        if not(len(ref) == 1 or len(alt) == 1):
            raise ValueError('Complex variant at %s' % pos)
        self.freq = altdepth / (refdepth + altdepth)
        self.span = len(ref) - len(alt)
        if not self.span:
            self.type = 'snv'
            self.repr = self.alt
            self.affected_pos = list(range(self.pos, self.pos+self.span+1))
        elif self.span > 0:
            self.type = 'del'
            self.repr = 'del'
            self.affected_pos = list(range(self.pos+1, self.pos+self.span+1))
        else:
            self.type = 'ins'
            self.repr = 'ins'
            self.affected_pos = list(range(self.pos, self.pos+self.span+1))

parser = argparse.ArgumentParser()

parser.add_argument('-p', '--pileup', help='Path to samtools mpileup output')
#parser.add_argument('-s', '--cns', help='Path to VarScan cns VCF')
parser.add_argument('--workflow_id')
parser.add_argument('-t', '--threshold', type=float, nargs='*', help='Threshold for consensus')
parser.add_argument('-i', '--iupac-thresh', type=float, help='Threshold for use of IUPAC code for split consensus between two alleles')
parser.add_argument('-o', '--output', help='Path to output TSV')
args = parser.parse_args()

#cns_df = pd.read_csv(args.cns, sep='\t', comment="#", header=None, names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'Sample'])

# Read variants from VarScan CNS file
"""
refseq = ''
exppos = 1
variants = []
for pos, ref, alt, data in zip(cns_df['POS'], cns_df['REF'], cns_df['ALT'], cns_df['Sample']):
    if exppos != int(pos):
        raise ValueError('Reference positions non-sequential')
    exppos += 1
    refseq += ref[0]
    if ',' in ref:
        raise ValueError('Multiple reference alleles')
    if ',' in alt:
        raise ValueError('Multiple alternate alleles')

    try:
        gt, gq, sdp, dp, rd, ad, freq, pval, rbq, abq, rdf, rdr, adf, adr = data.split(':')
    except ValueError: # Line will lack these fields if position is uncovered
        continue
    refdepth = int(rd)
    altdepth = int(ad)
#    if not refdepth + altdepth:
#        raise ValueError('Uncovered position %s' % pos)
    if alt != '.':
        var = Variant(pos, ref, alt, refdepth, altdepth)
        variants.append(var)

variant_lookup = {index: [var for var in variants if index in var.affected_pos] for index in range(1, len(refseq)+1)}
insertion_lookup = {index: [var for var in variants if (var.type == 'ins' and var.pos == index)] for index in range(1, len(refseq)+1)}
"""

if args.iupac_thresh < .334 or args.iupac_thresh > 0.50:
    raise ValueError('IUPAC code threshold must be between 0.334 and 0.50 inclusive')

alleles = ('A', 'C', 'G', 'T', 'del')

threshsort = sorted(args.threshold, reverse=True)
threshold_colnames = ['%s%% ID' % int(thresh*100) for thresh in threshsort]
#for thresh in args.threshold:
#    threshold_colnames += ['NGS Consensus %s' % thresh, 'Match %s' % thresh]

with open(args.output, 'w') as out:
    header_fields = ['Sample', 'Position', 'Reference', 'Count', 'Matches', 'Mismatches'] + list(alleles) + ['insert', 'consensus', '% Identity'] + threshold_colnames # + ['Variants', 'Insertions', 'Frequencies']
    if args.workflow_id:
        header_fields = ['workflow_id'] + header_fields
    out.write('\t'.join(header_fields) + '\n')

with open(args.pileup, 'r') as pileup:
    for linenum, line in enumerate(pileup):
        index = linenum + 1
        contig, pos, ref, count, basecalls, quality = line.rstrip().split('\t')
        count = int(count)
        insertion_count = len(list(re.finditer("\+[0-9]+", basecalls)))
        if args.workflow_id:
            outputs = [args.workflow_id]
        else:
            outputs = []
        outputs += [contig, pos, ref, count]
#        if not count:
#            frac = {allele: 0.0 for allele in alleles}
#        else:

        counts = {allele: 0.0 for allele in alleles}
            #frac = {}
        # Remove deletion starts from pileup (still present at covered loci)
        while True:
            dels = list(re.finditer("-[0-9]+", basecalls))
            if not dels:
                break
            match = dels[0]
            extent = int(match.group()[1:])
            basecalls = basecalls[:match.start()] + basecalls[match.end()+extent:]
        # Remove insertions from pileup
        while True:
            ins = list(re.finditer("\+[0-9]+", basecalls))
            if not ins:
                break
            match = ins[0]
            extent = int(match.group()[1:])
            basecalls = basecalls[:match.start()] + basecalls[match.end()+extent:]

        basecalls = re.sub("\^.", "", basecalls)
        basecalls = basecalls.replace("$", "").upper()
        if not count:
            basecalls = basecalls.replace("*", "")
        basecalls = basecalls.replace(".", ref)
        basecalls = basecalls.replace(",", ref)
        if count != len(basecalls):
        #    with open('debug.txt', 'w') as debug:
        #        debug.write(line)
            raise RuntimeError('Calls miscounted at line %s' % pos)
        for base in counts:
            counts[base] = basecalls.count(base)
        counts['del'] = basecalls.count('*')
        outputs += [counts[ref], count-counts[ref]]

        outputs += [str(counts[i]) for i in alleles]

        outputs.append(insertion_count)

        allele_fracs = {}
        for allele in alleles:
            divcount = max(count, 1)
            allele_fracs[allele] = counts[allele] / divcount

        consensus = None

        if allele_fracs['G'] >= args.iupac_thresh:
            if allele_fracs['C'] >= args.iupac_thresh:
                consensus = 'S'
            elif allele_fracs['A'] >= args.iupac_thresh:
                consensus = 'R'
        elif allele_fracs['T'] >= args.iupac_thresh:
            if allele_fracs['C'] >= args.iupac_thresh:
                consensus = 'Y'
            elif allele_fracs['A'] >= args.iupac_thresh:
                consensus = 'W'

        if not consensus:
            for allele in alleles:
                if allele_fracs[allele] >= max(args.threshold):
                    consensus = allele
                    break
            else:
                consensus = 'N'

        outputs.append(consensus)

        perc_ident = counts[ref] / max(count, 1)
        outputs.append("{:0.3f}".format(100*perc_ident))

        #strfrac = {k: "{:0.3f}".format(v) for k,v in frac.items()}

        # Get consensus at each threshold
        for thresh in threshsort:
#           for alt, freq in frac.items():
#                if freq > thresh:
#                    consensus = alt
#                    break
#            else:
#                consensus = 'N'
#            if consensus == ref:
#                match = 'Y'
#            else:
#                match = 'N'
            if perc_ident >= thresh:
                consensus_exists = 'Y'
            else:
                consensus_exists = 'N'

            outputs.append(consensus_exists)

        # Get called variants
        """
        var_at_pos = variant_lookup[index]
        observed_alts = [var.repr for var in var_at_pos] 
        outputs.append('-' + ','.join(observed_alts))
        ins_at_pos = insertion_lookup[index]
        inserted_sequences = [var.alt[1:] for var in ins_at_pos]
        insertion_freqs = ["{:0.3f}".format(var.freq) for var in ins_at_pos]
        outputs.append(','.join(inserted_sequences))
        outputs.append(','.join(insertion_freqs))
        """

        with open(args.output, 'a') as output:
            output.write('\t'.join([str(i) for i in outputs]) + '\n')
