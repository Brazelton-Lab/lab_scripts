#! /usr/bin/env python

# counts number of FASTA sequences and total bases for each file with given extension
# example usage:
# python fasta-summary.py .fa

import sys
extension = '*' + sys.argv[1]

from Bio import SeqIO
from Bio.SeqUtils import GC

l = []
import glob
for filename in glob.glob(extension):	
	l.append(filename)

for filename in l.sorted():
	seqs = 0
	bp = 0		
	for fasta in SeqIO.parse(filename,"fasta"):
		length = len(fasta.seq)
		seqs = seqs + 1
		bp = bp + int(length)
			
	print filename + ':',
	print str(seqs) + ' sequences and',
	print str(bp) + ' total bases'
