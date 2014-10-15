#! /usr/bin/env python

# extract list of seq names from raxml EPA classification file that have pendant lenghts > given value
# usage:
# python get-long-branches.py RAxML_classification.Hydrog.noOTUs.EPA1 0.1

import sys
cfilename = sys.argv[1]
cutoff = sys.argv[2]
outfilename = cfilename + '.longbranches.txt'
outfile = open(outfilename, 'a')

count = 0
cfile = open(cfilename)
for line in cfile:
	cols = line.split(' ')
	name = cols[0]
	length = cols[3].strip('\n')
	if float(length) > float(cutoff): 
		outfile.write(name+'\n')
		count = count + 1
	else: pass

print count,
print 'sequences with pendant length >',	
print cutoff	
	
cfile.close()
outfile.close()
