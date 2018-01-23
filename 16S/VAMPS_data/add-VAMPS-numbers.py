#! /usr/bin/env python

"""
reports total number of seqs in VAMPS fasta file that contains only unique seqs
the script just adds the numbers provided at the end of each fasta header
will process all files with given extension
usage:
python add-VAMPS-numbers.py .fa


Copyright:

    add-VAMPS-numbers  adds the numbers provided at the end of each fasta header

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
extension = '*' + sys.argv[1]

import glob
for filename in glob.glob(extension):	
	seqs = 0
	with open(filename) as file:
		for line in file:
			if line[0] == '>':
				words = line.split('|')
				num = int(words[-1].strip('\n'))
				seqs = seqs + num
	print filename,
	print seqs
