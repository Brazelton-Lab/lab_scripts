#! /usr/bin/env python

# extract sequences from FASTA file according to headers in provided file
# usage:
# python fasta-get.py file.fa names.txt

import sys
fastafilename = sys.argv[1]
namesfilename = sys.argv[2]
outfilename = fastafilename + '.select.fa'

l = []
with open(namesfilename) as namesfile:
	for name in namesfile: l.append(name.strip('\n'))
	
from Bio import SeqIO

with open(outfilename,'a') as outfile:
	for fasta in SeqIO.parse(fastafilename,'fasta'):
		if fasta.id in l: SeqIO.write(fasta,outfile,'fasta')	
