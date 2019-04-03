#! /usr/bin/env python

"""
reports summary of annotations in gff file

usage:
python prokka-counts.py <sample.gff>

Copyright:

    prokka-counts  reports summary counts from prokka gff file

    Copyright (C) 2019  William Brazelton

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
gff_file = sys.argv[1]

entries = 0
other = 0
none = 0
d = {}

with open(gff_file) as gff:
	for line in gff:
		if line[0] == '#': pass
		else: 
			entries = entries + 1
			if 'database' in line:
				fields = line.split(';')
				for field in fields:
					if 'database' in field: 
						db = field.split('=')
						db = db[1].strip('\n')
						if db in d: d[db] = d[db] + 1
						else: d[db] = 1
			elif 'no annotation' in line: none = none + 1
			else: other = other + 1

print str(entries) + '\t total entries in .gff'
for db in d:
	print str(d[db]) + '\t' + db
print str(other) + '\t entries with an annotation but no database specified'
print str(none) + '\t entries without an annotation'
