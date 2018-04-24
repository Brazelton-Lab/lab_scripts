#! /usr/bin/env python
# zero-out taxa for a user-provided specific sample 
# might be useful for removing taxa from some samples (but not others) prior to running sourcetracker, for example
# the script will find taxa in the provided count table that match the first column of a user-provided file
# assumes files are csv
# this version (distinct from zero-taxa.py) zeroes out all taxa except those provided by user

import sys
infilename = sys.argv[1]
listfilename = sys.argv[2]
sample = sys.argv[3]
outfilename = infilename.replace('.csv','.zeros.csv')

with open(listfilename) as listfile:
	l = []
	for row in listfile:
		row = row.split(',')
		i = row[0]
		l.append(i.strip())
	print l

with open(outfilename, 'w') as outfile:
	with open(infilename) as infile:
		for row in infile:			
			outfile.write(row)		# write the first line
			s = []
			cols = row.split(',')
			for col in cols:
				s.append(col)		# store list of sample names
			break
		
		if sample in s:
			pos = s.index(sample)
		else: 
			print 'error:',
			print sample, 
			print 'not found in',
			print infilename
			sys.exit()
			
		for row in infile:
			cols = row.split(',')
			j = cols[0]
			if j in l: outfile.write(row)
			else: 
				outfile.write(j)
				for col in cols[1:pos]: outfile.write(',' + col)
				outfile.write(',0')
				for col in cols[pos+1:]: outfile.write(',' + col)
			
