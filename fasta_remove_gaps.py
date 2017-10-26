#! /usr/bin/env python

# usage:
# python fasta_remove_gaps.py file.fa
# output: file.nogaps.fa

import sys
filename = sys.argv[1]

new = filename.split('.')
new = new[:-1]
outfilename = ''
for n in new: outfilename = outfilename + n + '.'
outfilename = outfilename + 'nogaps.fa'

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

with open(outfilename, 'w') as o:
	for fasta in SeqIO.parse(filename,"fasta"):
		f = Seq(str(fasta.seq), generic_dna)
		g = f.ungap("-")
		h = g.ungap(".")
		fasta.seq = h
		SeqIO.write(fasta,o,'fasta')
