#! /usr/bin/env python3

# searches lines of a gff file for user-provided text
# corresponding sequence name in FASTA file is written to new file
# also writes corresponding GFF file
# FASTA header must match the name in the first column of the gff file

import sys
import argparse

parser = argparse.ArgumentParser(description='searches gff file for keyword and writes corresponding FASTA sequence')
parser.add_argument('-g', '--gff', help='name of gff file to be searched')
parser.add_argument('-f', '--fasta', help='name of fasta file to be searched')
parser.add_argument('-t', '--text', help='text to search with')
parser.add_argument('-o', '--output', help='base name of fasta and gff output files')
args = parser.parse_args()

if args.output: 
	Foutfile = args.output + '.fa'
	Goutfile = args.output + '.gff'
if not args.output: 
	Foutfile = args.gff + 'select.fa'
	Goutfile = args.gff + 'select.gff'
	


l = []		# make list of sequence IDs whose gff entries contain the search term
with open(args.gff) as g, open(Goutfile, 'w') as Go:
	Go.write('##gff-version 3\n')
	for line in g:
		if args.text in line:
			cols = line.split('\t')
			ID = cols[0]
			l.append(ID)
			Go.write(line)
			
with open(args.fasta) as f, open(Foutfile, 'w') as Fo:
	status = 'no'
	for line in f:
		if line[0] == '>':
			header = line.split()
			header = header[0].strip('>')
			if header in l:
				#print(header)
				Fo.write('>')
				Fo.write(header)
				Fo.write('\n')
				status = 'yes'
			else: status = 'no'
				
		elif status == 'yes':
			Fo.write(line)
		
		else: pass
				
			
