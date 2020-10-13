#!/usr/bin/env python3

# searches lines of a gff file for user-provided text
# corresponding sequence name in FASTA file is written to new file
# FASTA header must match the name in the first column of the gff file

import sys
import argparse

parser = argparse.ArgumentParser(description='searches gff file for keyword and writes corresponding FASTA sequence')
parser.add_argument('-g', '--gff', help='name of gff file to be searched')
parser.add_argument('-f', '--fasta', help='name of fasta file to be searched')
parser.add_argument('-t', '--text', help='text to search with')
parser.add_argument('-o', '--output', help='name of fasta output file')
args = parser.parse_args()

if args.output: outfile = args.output
if not args.output: outfile = args.file + 'select.fa'

l = []		# make list of sequence IDs whose gff entries contain the search term
with open(args.gff) as g:
	for line in g:
		if args.text in line:
			cols = line.split('\t')
			ID = cols[0]
			l.append(ID)
			
with open(args.fasta) as f, open(outfile, 'w') as o:
	status = 'no'
	for line in f:
		
	
		if line[0] == '>':
			header = line.split()
			header = header[0].strip('>')
			if header in l:
				#print(header)
				o.write('>')
				o.write(header)
				o.write('\n')
				status = 'yes'
			else: status = 'no'
				
		elif status == 'yes':
			o.write(line)
		
		else: pass
				
			