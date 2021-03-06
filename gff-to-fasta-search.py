#! /usr/bin/env python3

# searches lines of a gff file for user-provided text
# corresponding sequence name in FASTA file is written to new file
# also writes corresponding GFF file
# FASTA header must match the name in the ID field
# if FASTA headers contain two underscores, e.g. c_000000000001_1, these are assumed to represent protein sequences (e.g. .faa file) and the headers must match the ID field exactly
# if FASTA headers contain one undersacore, e.g. c_000000000001, these are assumed to represent contig DNA sequences (e.g. .fa file) and the headers must match the ID field, with the second underscore and following characters omitted

import sys
import argparse

parser = argparse.ArgumentParser(description='searches gff file for keyword and writes corresponding FASTA sequence')
parser.add_argument('-g', '--gff', help='name of gff file to be searched')
parser.add_argument('-f', '--fasta', help='name of fasta file to be searched; can be either contig DNA or individual protein sequences')
parser.add_argument('-t', '--text', help='text to search with')
parser.add_argument('-o', '--output', help='base name of fasta and gff output files')
parser.add_argument('-m', '--mode', help='mode: options are "contig" (FASTA header must match first column of gff line; i.e. contig name) or "gene" (FASTA header must match ID field; i.e. ORF name)')
args = parser.parse_args()

if args.output:
	if args.mode == 'contig': 
		Foutfile = args.output + '.fa'
		Goutfile = args.output + '.gff'
	elif args.mode == 'gene':
		Foutfile = args.output + '.faa'
		Goutfile = args.output + '.gff'
elif not args.output: 
	if args.mode == 'contig':
		Foutfile = args.gff + 'select.fa'
		Goutfile = args.gff + 'select.gff'
	elif args.mode == 'gene':
		Foutfile = args.gff + 'select.faa'
		Goutfile = args.gff + 'select.gff'

l = []	# make list of sequence IDs whose gff entries contain the search term
with open(args.gff) as g, open(Goutfile, 'w') as Go:
	Go.write('##gff-version 3\n')
	for line in g:
		if args.text in line:
			if args.mode == "gene":
				cols = line.split('ID=')
				ID = cols[1].split(';')
				l.append(ID[0])
				Go.write(line)
			elif args.mode == "contig":
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
				print(header)
				Fo.write('>')
				Fo.write(header)
				Fo.write('\n')
				status = 'yes'
			else: status = 'no'
		elif status == 'yes':
			Fo.write(line)
		
		else: pass
				
			
