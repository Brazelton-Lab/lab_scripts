# sort FASTA entries according to their abundance provided in .names file
# typically used after the unique.seqs command in mothur
# usage:
# python sort-fasta-according-to-names.py filename.unique.fa filename.names

import sys
fasta_filename = sys.argv[1]
names_filename = sys.argv[2]
outfilename = fasta_filename + '.sorted.fa'

l = []
with open(names_filename) as names_file:
	for line in names_file:
		name = line.split('\t')
		name = name[0]
		length = int(len(line.split(',')))
		l.append([length,name])

l.sort()
l.reverse()

from Bio import SeqIO
with open(outfilename,'a') as outfile:
	for i in l:
		for fasta in SeqIO.parse(fasta_filename,"fasta"):
			name = fasta.id
			if name == i[1]: 
				SeqIO.write(fasta,outfile,"fasta")
				break	
			else: pass
print len(l),
print 'unique sequences sorted from most to least abundant'	
