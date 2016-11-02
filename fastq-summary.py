#! /usr/bin/env python

# counts number of FASTA sequences and total bases for each file with given extension
# example usage:
# python fasta-summary.py .fa

"""
Copyright:

    fastq-summary.py Report number of sequences and bases in FASTQ file
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

from Bio import SeqIO
from Bio.SeqUtils import GC

import glob
for filename in glob.glob(extension):
	seqs = 0
	bp = 0
	for fastq in SeqIO.parse(filename,"fastq"):
		length = len(fastq.seq)
		seqs = seqs + 1
		bp = bp + int(length)

	print filename + ':',
	print str(seqs) + ' sequences and',
	print str(bp) + ' total bases'
