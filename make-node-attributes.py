#! /usr/bin/env python

"""
requires two input files: 
1. a csv file of the abundance of each sequence in each sample, with sequences as rows. You probably want transformed counts, i.e. proportions, not raw counts 
2. a mothur taxonomy file, preferably after reformatting it with taxonomy_edit.py

usage:
make-node-attributes.py proportions.csv filename.taxonomy.renamed.txt outfilename.txt
outfilename is optional. default is node-attributes.txt

Copyright:

    make-node-attributes  writes both the sample name and the abundance for that seq name separated by tabs

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
taxfilename = sys.argv[1]
propfilename = sys.argv[2]
try: outfilename = sys.argv[3]
except: 
	outfilename = 'node-attributes.txt'
	with open(outfilename, 'w') as outfile: 	# if node-attributes.txt already exists, erase it
		outfile.write('Sequence\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tStrain\tMaxSample\tMaxAbundnace\n')


with open(propfilename) as propfile:
	for line in propfile:
		line = line.strip('\n')
		line = line.replace('"','')
		headers = line.split(',')
		headers = headers[1:]	# need list of sample names for later
		break	# breaks after first line because we only need the sample names

	d = {}
	for line in propfile:		# did not re-open the file so that the for: loop skips the first line
		line = line.strip('\n')
		cols = line.split(',')
		seqname = cols[0].replace('"','')
		numbers = []
		for i in cols[1:]: numbers.append(float(i))
		maxcolumn = numbers.index(max(numbers))
		d[seqname] = [headers[maxcolumn],max(numbers)]	# assigns a list of sample name and abundance to the dictionary
        
with open(taxfilename) as taxfile:
	yes = 0
	no = 0
	for line in taxfile:
		line = line.strip('\n')
		cols = line.split('\t')
		taxname = cols[0]
		if taxname in d:
			yes = yes + 1
			with open(outfilename,'a') as outfile:
				line = line.rstrip('\t')		# removes last tab to avoid two tabs when dictionary items are written to file
				outfile.write(line)
				for i in d[taxname]:			# looks up maxcolum assigned to that seq name in dictionary 
					outfile.write('\t'+str(i))		# writes both the sample name and the abundance for that seq name separated by tabs
				outfile.write('\n')       
		else: no = no +1

print 'found abundances for',
print yes,
print 'sequences'
print 'did not find abundances for',
print no,
print 'sequences'				
