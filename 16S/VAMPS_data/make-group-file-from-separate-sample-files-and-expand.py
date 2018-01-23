#! /usr/bin/env python
"""

make group file where the sample name is the fasta filename, not in the fasta header
in other words, takes a folder full of .fasta files and makes one group file
also expands the file names as in fasta-expander-vamps-2014June.py
specify the extension in the command
script must be run from same folder as the files, or else specify the path in the script below

usage example 1:
python make-group-file-from-separate-sample-files.py .fa

usage example 2:
python make-group-file-from-separate-sample-files.py .unique



Copyright:

    make-group-file-from-separate-sample-files-and-expand  make group file where the sample name is the fasta filename, not in the fasta header

    Copyright (C) 2016  William Brazelton <comma-separated list of authors>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

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
		abundance = fasta.description.split('frequency:')
		abundance = abundance[1]
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
