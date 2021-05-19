#! /usr/bin/env python3

# only keep rows that contain text provided by user
# example
# python filter-table.py -i LC2018_all_fluids.IDabunds.mg.mt.tpm.txt.KEGG.tax.tsv -k KO-keep.txt -c 1 -o LC2018_all_fluids.IDabunds.mg.mt.tpm.txt.KEGG-keep.tax.tsv

import sys
import argparse

parser = argparse.ArgumentParser(description='only keep rows that contain text provided by user')
parser.add_argument('-i', '--input', help='input file')
parser.add_argument('-k', '--keep', help='list of items to keep if they appear anywhere in each line of input file')
parser.add_argument('-c', '--column', help='number of column to search for items to keep, where the first column = 0')
parser.add_argument('-o', '--output', help='name of new file')
args = parser.parse_args()

# create list of items to keep
L = []
with open(args.keep) as k:
	for line in k:
		L.append(line.strip('\n'))

# iterate through input and write output
with open(args.input) as i, open(args.output, 'w') as o:
	for line in i:
		cols = line.split('\t')
		if len(cols) < 2: cols = line.split(',')
		col = cols[int(args.column)]
		if col in L: o.write(line)
		else: pass
