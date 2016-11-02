#!/usr/bin/env python
# given output from pathways2contigs.py, filter table according to coverage or pathway
# usage:
# python filter-pathways-table.py file.IDs.contigs.tsv file.IDs.contigs.filter.tsv 0.5 all
# change 0.5 to the desired coverage threshold (inclusive, pathway with coverage of exactly 0.5 will be written to output file)
# change all to a specific pathway - only one allowed, but it can be general or specific
# i.e. this is valid input: Transporters;ABC_transporters,_prokaryotic_type;ABC-2_type_and_other_transporters;Antibiotic_ABC_transporter_[mD:m00248]
# and this is valid input: Transporters

"""
Copyright:

    filter-pathways-table.py Filter pathways table by coverage
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
newfilename = sys.argv[2]
coverage = sys.argv[3]
pathway = sys.argv[4]

if "'" in pathway: pass
else: print 'WARNING: pathway search is not accurate unless search term is enclosed in quotes'

l = 0
yes = 0
no = 0
with open(newfilename,'w') as newfile:
	with open(filename) as file:
		for line in file:
			l = l + 1
			cols = line.split('\t')
			if float(cols[1]) < float(coverage): no = no + 1
			else:
				if pathway == 'all':
					newfile.write(line)
					yes = yes + 1
				elif pathway in line:
					newfile.write(line)
					yes = yes + 1
				else: no = no + 1

print str(l) + ' lines in table'
print str(yes) + ' lines written to ' + newfilename
print str(no) + ' lines did not match criteria'

