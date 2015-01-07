#! /usr/bin/env python
# partition FASTA file into multiple files according to header

import sys
infilename = sys.argv[1]

from Bio import SeqIO
for seq in SeqIO.parse(infilename,"fasta"):
	words = seq.id.split('|')
	sample = words[2]
	outfilename = sample + '.fa'
	with open(outfilename, 'a') as outfile:
		SeqIO.write(seq,outfile,"fasta")