#! /usr/bin/env python

# adds "X. to the front and ." to the end of each name in the first column if it begins with Otu
# used to format network attributes file for Cytoscape so that Otu names match
# assumes comma separated values and unix linebreaks
# usage:
# python format-for-cytoscape.py filename.csv

"""
Copyright:

    format-for-cytoscape.py Modify OTU file for use with Cytoscape
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
outfilename = filename.replace('.csv','.for-cys.csv')

status = 'no'
with open(filename) as file:
	for line in file:
		columns = line.split(',')
		if 'Otu' in columns[0]:
			status = 'yes'
			newname = '"X.' + columns[0] + '."'
			with open(outfilename,'a') as outfile:
				outfile.write(newname)
				for col in columns[1:]:
					outfile.write(',' + col)
		else:
			with open(outfilename,'a') as outfile:
				outfile.write(line)

if status == 'no': print 'Nothing happened. Perhaps your input file needs to be converted to Unix linebreaks? Delete the output file and try again.'
