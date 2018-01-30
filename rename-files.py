#! /usr/bin/env python
"""
rename files according to provided table, which should be tab-delimited
directory of files to be renamed must be specified. if files are in same directory, enter: ./

usage:
python rename-files.py table-with-conversions.txt directory

Copyright:

    rename-files  rename files according to provided table

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
    along with this program.  If not, see <http://www.gnu.org/licenses/>.﻿
"""

import sys
infilename = sys.argv[1]
dirname = sys.argv[2]

D = {}	# make dictionary of names and ids
with open(infilename) as infile:
	for line in infile:
		cols = line.split('\t')
		name = cols[0]
		id = cols[1].strip('\n')
		D[id] = name
		
import os
from subprocess import call
for root, dir, files in os.walk(dirname):
	for file in files:		
		print file
		filename_id = file.split('.')
		filename_id = filename_id[0]
		if filename_id in D: 
			oldname = os.path.join(root,file)
			newname = oldname.replace(filename_id,D[filename_id])
			call(['mv',oldname,newname])
		else: pass	
		
