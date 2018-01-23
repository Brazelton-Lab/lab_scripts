#! /usr/bin/env python

"""

Copyright:

    transpose-table  moves data in the rows of a file to the columns separated by tabs

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
outfilename = infilename + '.transposed'

with open(infilename) as infile:
	matrix = []
	for line in infile:
		items = line.split('\t')
		L = []
		for item in items:
			item = item.strip('\n')
			L.append(item)
		matrix.append(L)

#print matrix	
transposed = map(None,*matrix)
#print transposed

with open(outfilename, 'a') as outfile:
	for line in transposed:
		for item in line:
			if item == None: pass
			else: 
				outfile.write(item)
				outfile.write('\t')
		outfile.write('\n')	
