#! /usr/bin/env python

"""


Copyright:

    fasta-expander-vamps-2016  expand VAMPS fasta file according to abundance indicate at end of header

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
filename = sys.argv[1]

outfilename = filename.replace('.fa','.expanded.fa')
outfile = open(outfilename,'a')

from Bio import SeqIO

for fasta in SeqIO.parse(filename,"fasta"):
	d = fasta.description
	root = fasta.id
	abundance = d.split('|')
	abundance = abundance[-1]
	if int(abundance) == 1:	SeqIO.write(fasta,outfile,"fasta")	
	else:
		count = 1
		SeqIO.write(fasta,outfile,"fasta")
		while 1:
			count = count + 1
			fasta.id = root + 'r' + str(count)
			#print fasta.id
			SeqIO.write(fasta,outfile,"fasta")
			if int(count) == int(abundance): break
		
outfile.close()		
