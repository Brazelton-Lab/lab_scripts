#! /usr/bin/env python
# edit user-provided ESOM class file with new assignments in user-provided file
# each line of user-provided file of new assignments should contain a data point number and a class number, separated by tabs
# usage:
# python edit-esom-class-file.py esom.cls new-assignments.tsv new-class-filename.cls

"""
Copyright:

    edit-esom-class-file.py Append user data to ESOM class file
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
cls_file = sys.argv[1]
user_file = sys.argv[2]
new_file = sys.argv[3]

# create dictionary of user-provided new assignments:
d = {}
with open(user_file) as user:
	for line in user:
		cols = line.split('\t')
		data_point = cols[0].strip()
		cls_number = cols[1].strip()
		d[data_point] = cls_number.strip('\n')

# iterate through class file, writing new class file with new assignments:
with open(new_file,'w') as new:
	with open(cls_file) as cls:
		for line in cls:
			if line[0] == '%': new.write(line)
			else:
				cols = line.split('\t')
				if cols[0] in d: new.write(str(cols[0]) + '\t' + str(d[cols[0]]) + '\n')
				else: new.write(line)
print 'WARNING: if you introduced new classes to this .cls file, you need to manually add them to the header of this new .cls file'





