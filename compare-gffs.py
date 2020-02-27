#! /usr/bin/env python

# scans all .gff files in current directory for the presence/absence of IDs provided by the user
# generates two output files, one with gene IDs as columns and one with gene IDs as rows
# searches for the user-provided gene ID in all text in the .gff file. not field-specific.

import sys
import argparse

parser = argparse.ArgumentParser(description='scan .gff files for provided list of gene IDs')

parser.add_argument('-t','--tableFile', help='table (or list) of IDs to search for', required=True)
parser.add_argument('-o','--outFile', help='output table file', required=True)
args = parser.parse_args()

l = []
with open(args.tableFile) as table:
	for line in table: 
		line = line.replace('\t',',')
		cols = line.split(',')
		id = cols[0].strip('\n')
		l.append(id)

with open(args.outFile, 'w') as outfile:
	# write header row of new file
	outfile.write('IDs')
	for i in l:
		outfile.write(',')
		outfile.write(i)
	outfile.write('\n')
	
	import glob
	for filename in glob.glob('*.gff*'):
		with open(filename) as gff:
			contents = gff.read()
			outfile.write(filename)
			for i in l:
				outfile.write(',')
				if i in contents: outfile.write(i)
				else: pass
			outfile.write('\n')

#print('finished writing',args.outFile)

# output has gene IDs as columns, so now open that output file, transpose it, and save as new file
import pandas
trfilename = args.outFile + '.tr.csv'
t1 = pandas.read_csv(args.outFile, delimiter=',', low_memory=False)
t2 = t1.transpose()
t2.to_csv(trfilename, sep=',', header=False)

#print('finished writing',trfilename)
