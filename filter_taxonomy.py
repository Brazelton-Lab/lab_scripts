#!/usr/bin/env python

# select lines from .taxonomy file that match the OTU name in a count_table file
# usage:
# filter_taxonomy.py yourfile.taxonomy yourfile.count_table

import sys
taxfilename = sys.argv[1]
countfilename = sys.argv[2]
newfilename = taxfilename + '.filtered'

l = []
with open(newfilename,'w') as newfile:
	with open(countfilename) as cfile:
		for line in cfile:
			cols = line.split('\t')
			l.append(cols[0])
	with open(taxfilename) as tfile:
		for line in tfile:
			cols = line.split('\t')
			if cols[0] in l: newfile.write(line)			
			else: pass
