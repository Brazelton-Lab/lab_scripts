# reports total number of seqs in VAMPS fasta file that contains only unique seqs
# the script just adds the numbers provided at the end of each fasta header
# will process all files with given extension
# usage:
# python add-VAMPS-numbers.py .fa

import sys
extension = '*' + sys.argv[1]

import glob
for filename in glob.glob(extension):	
	seqs = 0
	with open(filename) as file:
		for line in file:
			if line[0] == '>':
				words = line.split('|')
				num = int(words[-1].strip('\n'))
				seqs = seqs + num
	print filename,
	print seqs
