#! /usr/bin/env python

"""
parses output from my combine_matrix_pqvalues.py script
exports rows that have p and q values both greater than 0.05

usage:
python parse-Qvalues.py your_cors_and_pqvalues.txt

Copyright:

    parse-Qvalues  parses output from the combine_matrix_pqvalues.py script

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
filename = 'your_cors_and_pqvalues.txt'
f = open(filename)

outfilename = filename.replace('.txt','.05.txt')
outfile = open(outfilename, 'a')

for line in f:
	columns = line.split('\t')
	if columns[0] == 'var1': outfile.write(line)	# always writes first line
	elif len(columns) < 2: pass	# skips last line
	else:	
		if float(columns[4]) < 0.05:
			if float(columns[5]) < 0.05:
				outfile.write(line)
 
f.close()
outfile.close() 
 
