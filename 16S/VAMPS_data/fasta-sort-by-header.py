#! /usr/bin/env python

"""


Copyright:

    fasta-sort-by-header  partition FASTA file into multiple files according to header

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
infilename = sys.argv[1]

from Bio import SeqIO
for seq in SeqIO.parse(infilename,"fasta"):
	words = seq.id.split('|')
	sample = words[2]
	outfilename = sample + '.fa'
	with open(outfilename, 'a') as outfile:
		SeqIO.write(seq,outfile,"fasta")
