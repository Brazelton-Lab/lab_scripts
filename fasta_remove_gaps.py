#! /usr/bin/env python


"""

usage:
python fasta_remove_gaps.py file.fa
output: file.nogaps.fa

Copyright:

    fasta_remove_gaps  Removes gaps from a .fasta file

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
filename = sys.argv[1]

new = filename.split('.')
new = new[:-1]
outfilename = ''
for n in new: outfilename = outfilename + n + '.'
outfilename = outfilename + 'nogaps.fa'

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

with open(outfilename, 'w') as o:
	for fasta in SeqIO.parse(filename,"fasta"):
		f = Seq(str(fasta.seq), generic_dna)
		g = f.ungap("-")
		h = g.ungap(".")
		fasta.seq = h
		SeqIO.write(fasta,o,'fasta')
