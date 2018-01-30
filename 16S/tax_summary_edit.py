#! /usr/bin/env python

"""
edits mothur taxonomy summary file
transfers last name that is not "unclassified" or "uncultured" to "unclassified" or "uncultured" assignment
make sure that the file has default sorting (by rankID)

Copyright:

    tax_summary_edit  edits mothur taxonomy summary file

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
outfilename = infilename + '.renamed.txt'
outfile = open(outfilename,'a')

infile = open(infilename)
for line in infile:
	if "unclassified" in line:
		columns = line.split('\t')
		tax = columns[2]
		newtax = tax + ' ' + lasttax
		outfile.write(columns[0])
		outfile.write('\t')
		outfile.write(columns[1])
		outfile.write('\t')
		outfile.write(newtax)
		for tab in columns[3:]:
			outfile.write('\t')
			outfile.write(tab)
	elif "uncultured" in line:
		columns = line.split('\t')
		tax = columns[2]
		newtax = tax + ' ' + lasttax
		outfile.write(columns[0])
		outfile.write('\t')
		outfile.write(columns[1])
		outfile.write('\t')
		outfile.write(newtax)
		for tab in columns[3:]:
			outfile.write('\t')
			outfile.write(tab)
	else: 
		outfile.write(line)
		columns = line.split('\t')
		lasttax = columns[2]
infile.close()
outfile.close()		
