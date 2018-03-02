#! /usr/bin/env python

"""
THIS SCRIPT NOT NECESSARY - JUST CLICK ON THE COLUMN IN CYTOSCAPE IMPORT WINDOW THAT YOU WANT TO BE THE EDGE ATTRIBUTES

from file with variable names and correlation coefficient, export Edge Attributes file for importing with Cytoscape
infile must end with ".txt"

usage: python write_edge_attributes.py filename

Copyright:

    write_edge_attributes  export Edge Attributes file for importing with Cytoscape

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
infilename = sys.argv[1]
outfilename = infilename.replace('.txt','_edge_attr.txt')
infile = open(infilename)
outfile = open(outfilename,'a')
outfile.write('Continuous edge Attr1\n')

for line in infile.xreadlines():
	if line[:4] == "var1": pass
	elif len(line) < 2: pass
	else:
		columns = line.split('\t')
		outfile.write(columns[0])
		if float(columns[2]) > 0: outfile.write(' (pos) ')
		elif float(columns[2]) < 0: outfile.write(' (neg) ')
		elif float(columns[2]) == 0: outfile.write(' (neutral) ')
		else: outfile.write(' (error) ')
		outfile.write(columns[1])
		outfile.write(' = ')
		outfile.write(columns[2])
		outfile.write('\n')
infile.close()
outfile.close()
