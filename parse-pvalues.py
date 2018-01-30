#! /usr/bin/env python

"""
parses p values output from rcor.test in R

usage:
python parse-pvalues.py yourfilename.rcortest.pvalues

Copyright:

    parse-pvalues  parses p values from rcor.test in R

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
filename = 'myfilename.rcortest.pvalues.csv'
f = open(filename)

outfilename = filename.replace('csv','txt')
outfile = open(outfilename, 'a')

for line in f:
	columns = line.split(',')
	if columns[0] == '""': pass	# skips first line
	elif len(columns) < 2: pass	# skips last line
	else:	
		pvalue = columns[3]
		if pvalue == '"NA\n"': print 'Your input data probably includes columns containing only zeroes. You should delete these and repeat the R commands.'
		else: outfile.write(pvalue)
 
f.close()
outfile.close() 
 
