#! /usr/bin/env python

"""
converts OTU file from QIIME into .list file for mothur
incorporates sequence names from provided .names file

usage: python convert_qiime_otus_to_mothur_list.py otus.txt VAMPS_JGIbtrim30.unique.good.filter.names 31776
user must provide the following input manually:


Copyright:

    convert_qiime_otus_to_mothur_list_with_names  converts OTU file from QIIME into .list file for mothur

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

dist = '0.03'

import sys
infilename = sys.argv[1]
namesfilename = sys.argv[2]
otus = sys.argv[3]
outfilename = infilename + '.list'
outfile = open(outfilename, 'a')
outfile.write(dist+'\t'+otus+'\t')

namesfile = open(namesfilename)
N=[]
for nline in namesfile:
	nline = nline.strip('\n')
	ntabs = nline.split('\t')
	names = ntabs[1]
	#names = names.split(',')
	N.append(names)	

infile = open(infilename)
for line in infile:
	tabs = line.split('\t')
	seq = tabs[1].strip('\n')
	#print seq
	found = 'no'
	for names in N:
		#print names
		if seq in names: 
			outfile.write(names)
			found = 'yes'
	if found == 'no': 
		print seq,
		print 'not found in names file'		
	if len(tabs) >2:
		#print tabs
		found = 'no'
		for tab in tabs[2:]:
			outfile.write(',')
			seq=tab.strip('\n')
			for names in N:
				if seq in names: 
					outfile.write(names)
					found = 'yes'
					#print names
			if found == 'no': 
				print seq,
				print 'not found in names file'			
	outfile.write('\t')
infile.close()
outfile.close()
namesfile.close()
