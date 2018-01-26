#! /usr/bin/env python

"""
extract list of seq names from raxml EPA classification file that have pendant lenghts > given value

usage:
python get-long-branches.py RAxML_classification.Hydrog.noOTUs.EPA1 0.1

Copyright:

    get-long-branches  extract list of seq names from raxml EPA classification file that have pendant lenghts > given value

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
cfilename = sys.argv[1]
cutoff = sys.argv[2]
outfilename = cfilename + '.longbranches.txt'
outfile = open(outfilename, 'a')

count = 0
cfile = open(cfilename)
for line in cfile:
	cols = line.split(' ')
	name = cols[0]
	length = cols[3].strip('\n')
	if float(length) > float(cutoff): 
		outfile.write(name+'\n')
		count = count + 1
	else: pass

print count,
print 'sequences with pendant length >',	
print cutoff	
	
cfile.close()
outfile.close()
