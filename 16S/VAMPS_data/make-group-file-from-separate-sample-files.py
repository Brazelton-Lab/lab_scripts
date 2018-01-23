#! /usr/bin/env python

"""
make group file where the sample name is the fasta filename, not in the fasta header
in other words, takes a folder full of .fasta files and makes one group file
must specify the extension in the command
script must be run from same folder that contains the files, or else specify the path in the script below

usage example 1:
python make-group-file-from-separate-sample-files.py .fa

usage example 2:
python make-group-file-from-separate-sample-files.py .unique


Copyright:

    make-group-file-from-separate-sample-files  make group file where the sample name is the fasta filename, not in the fasta header

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

outfilename = 'all.group'

path = r'./'
for infilename in glob.glob(os.path.join(path, '*' + extension)):	
	for fasta in SeqIO.parse(infilename,"fasta"):		
		with open(outfilename,'a') as outfile: 
			outfile.write(fasta.id + '\t')	
			
			sample = infilename.replace(extension,'')
			outfile.write(sample.strip('./') + '\n')
