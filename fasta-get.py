#! /usr/bin/env python


"""
extract sequences from FASTA file according to headers in provided file
usage:
python fasta-get.py file.fa names.txt

Copyright:

    fasta-get.py Extract FASTA sequences based on entry headers
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
fastafilename = sys.argv[1]
namesfilename = sys.argv[2]
outfilename = fastafilename + '.select.fa'

l = []
with open(namesfilename) as namesfile:
	for name in namesfile: l.append(name.strip('\n'))

from Bio import SeqIO

with open(outfilename,'a') as outfile:
	for fasta in SeqIO.parse(fastafilename,'fasta'):
		if fasta.id in l: SeqIO.write(fasta,outfile,'fasta')
