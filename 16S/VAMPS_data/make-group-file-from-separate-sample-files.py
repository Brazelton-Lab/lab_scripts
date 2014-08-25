# make group file where the sample name is the fasta filename, not in the fasta header
# in other words, takes a folder full of .fasta files and makes one group file
# specify the extension in the command
# script must be run from same folder that contains the files, or else specify the path in the script below

# usage example 1:
# python make-group-file-from-separate-sample-files.py .fa

# usage example 2:
# python make-group-file-from-separate-sample-files.py .unique

import sys
import os, glob
from Bio import SeqIO
extension = sys.argv[1]

outfilename = 'all.group'

path = r'./'
for infilename in glob.glob(os.path.join(path, '*' + extension)):	
	for fasta in SeqIO.parse(infilename,"fasta"):		
		with open(outfilename,'a') as outfile: 
			outfile.write(fasta.id + '\t')	
			
			sample = infilename.replace(extension,'')
			outfile.write(sample.strip('./') + '\n')
