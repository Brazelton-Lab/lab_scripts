#! /usr/bin/env python
# make group file where the sample name is the fasta filename, not in the fasta header
# in other words, takes a folder full of .fasta files and makes one group file
# also expands the file names as in fasta-expander-vamps-2013.py
# specify the extension in the command
# script must be run from same folder as the files, or else specify the path in the script below

# usage example 1:
# python make-group-file-from-separate-sample-files.py .fa

# usage example 2:
# python make-group-file-from-separate-sample-files.py .unique

import sys
import os, glob
from Bio import SeqIO
extension = sys.argv[1]

outfile = open('all.group','a')

path = r'./'
for infilename in glob.glob(os.path.join(path, '*' + extension)):	
	print infilename
	for fasta in SeqIO.parse(infilename,"fasta"):		
		sample = infilename.replace(extension,'')
		root = fasta.id
		abundance = fasta.description.split('|')
		abundance = abundance[-1]
		if int(abundance) == 1: outfile.write(fasta.id + '\t' + sample.strip('./') + '\n')	
		else:
			outfile.write(fasta.id + '\t' + sample.strip('./') + '\n') # first, original sequence before write the replicates
			count = 1
			while 1:
				count = count + 1
				fasta.id = root + 'r' + str(count)
				outfile.write(fasta.id + '\t' + sample.strip('./') + '\n')					
				if int(count) == int(abundance): break
outfile.close()			
