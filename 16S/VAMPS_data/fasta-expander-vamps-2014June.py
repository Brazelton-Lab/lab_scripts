#! /usr/bin/env python
# expand VAMPS fasta file according to abundance indicate at end of header

import sys
filename = sys.argv[1]

outfilename = filename + '.expanded.fa'
outfile = open(outfilename,'a')

from Bio import SeqIO

for fasta in SeqIO.parse(filename,"fasta"):
	root = fasta.id
	abundance = fasta.description.split('frequency:')
	abundance = abundance[1]
	if int(abundance) == 1:	SeqIO.write(fasta,outfile,"fasta")	
	else:
		count = 1
		SeqIO.write(fasta,outfile,"fasta")
		while 1:
			count = count + 1
			fasta.id = root + 'r' + str(count)
			#print fasta.id
			SeqIO.write(fasta,outfile,"fasta")
			if int(count) == int(abundance): break
		
outfile.close()		
