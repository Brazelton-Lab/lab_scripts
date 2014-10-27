# output is a new file with length and GC for each contig provided in input file
# will process any file with given extension
# example usage:
# python fasta-length-GC.py .fa

import sys
extension = '*' + sys.argv[1]

import glob
for filename in glob.glob(extension):
	outfilename = filename.replace(sys.argv[1],'.length-GC.txt')
	print outfilename
	
	from Bio import SeqIO
	from Bio.SeqUtils import GC

	for fasta in SeqIO.parse(filename,"fasta"):
		name = fasta.id
		length = len(fasta.seq)
		gc = GC(fasta.seq)
		
		with open(outfilename,'a') as outfile: 
			outfile.write(name + '\t' + str(length) + '\t' + str(gc) + '\n')