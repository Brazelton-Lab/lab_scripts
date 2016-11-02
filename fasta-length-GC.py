#! /usr/bin/env python

# output is a new file with length and GC for each contig provided in input file
# will process any file with given extension
# example usage:
# python fasta-length-GC.py .fa

"""
Copyright:

    fasta-length-GC.py Get length and GC of sequences per entry from FASTA file
    Copyright (C) 2016  William Brazelton

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
