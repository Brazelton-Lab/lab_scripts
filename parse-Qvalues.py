#! /usr/bin/env python

# parses output from my combine_matrix_pqvalues.py script
# exports rows that have p and q values both greater than 0.05
# usage:
# python parse-Qvalues.py your_cors_and_pqvalues.txt

import sys
filename = 'your_cors_and_pqvalues.txt'
f = open(filename)

outfilename = filename.replace('.txt','.05.txt')
outfile = open(outfilename, 'a')

for line in f:
	columns = line.split('\t')
	if columns[0] == 'var1': outfile.write(line)	# always writes first line
	elif len(columns) < 2: pass	# skips last line
	else:	
		if float(columns[4]) < 0.05:
			if float(columns[5]) < 0.05:
				outfile.write(line)
 
f.close()
outfile.close() 
 
