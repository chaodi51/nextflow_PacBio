#!/usr/bin/env python
  
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input')
parser.add_argument('--output')
args = parser.parse_args()

counts = {} 
with open(args.input, 'r') as infile:
    counter = 0
    for line in infile:
        if not counter:
            with open(args.output, 'w') as output:
                output.write(line)
            counter += 1
            continue
        name, pos, count = line.rstrip().split('\t')
        counts[int(pos)] = int(count)

for index in range(1, max(counts.keys()) + 1):
    if index not in counts:
        counts[index] = 0
    with open(args.output, 'a') as output:
        output.write('\t'.join([name, str(index), str(counts[index])]) + '\n')


