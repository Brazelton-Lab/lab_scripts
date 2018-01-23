#! /usr/bin/env python

"""
renames FASTA file according to .groups file
use when sending FASTA file from mothur to Qiime
usage:
python FASTA-rename-with-groups.py filename.fa filename.group


Copyright:

    FASTA-rename-with-groups  renames FASTA file according to .groups file

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
fastafilename = sys.argv[1]
groupfilename = sys.argv[2]

outfilename = ''
splitname = fastafilename.split('.')
for word in splitname[:-1]:
	outfilename = outfilename + word + '.'
outfilename = outfilename + 'groups.fa'

D = {}
with open(groupfilename) as groupfile:
	for line in groupfile:
		columns = line.split('\t')
		D[columns[0]] = columns[1].strip('\n')

count = 0
with open(fastafilename) as fasta:
	for line in fasta:
		if '>' in line:
			count = count + 1
			name = line.strip('\n')
			name = name.replace('>','')
			try: 
				groupname = D[name]
				with open(outfilename,'a') as outfile: outfile.write('>' + groupname + '_' + str(count) + ' ' + name + '\n')
			except: print name + ' not found'
		else: 
			with open(outfilename,'a') as outfile: outfile.write(line)	

	
