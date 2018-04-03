#! /usr/bin/env python
# export rows from a table where the first column matches the first column of a user-provided file
# for example: cross-reference sequences in count tables
# assumes files are csv

import sys
infilename = sys.argv[1]
listfilename = sys.argv[2]
outfilename = infilename.replace('.csv','.select.csv')

with open(listfilename) as listfile:
	l = []
	for row in listfile:
		row = row.split(',')
		i = row[0]
		l.append(i.strip('\n'))

with open(outfilename, 'w') as outfile:
	with open(infilename) as infile:
		for row in infile:			# write the first line
			outfile.write(row)
			break
		for row in infile:
			cols = row.split(',')
			j = cols[0]
			if j in l: outfile.write(row)
			else: pass
