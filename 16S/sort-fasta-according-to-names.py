#! /usr/bin/env python

"""
sort FASTA entries according to their abundance provided in .names file
typically used after the unique.seqs command in mothur

usage:
python sort-fasta-according-to-names.py filename.unique.fa filename.names

Copyright:

    sort-fasta-according-to-names  sort FASTA entries according to their abundance provided in .names file

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
