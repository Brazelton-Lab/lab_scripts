# rename files according to provided table, which should be tab-delimited
# directory of files to be renamed must be specified. if files are in same directory, enter: ./
# usage:
# python rename-files.py table-with-conversions.txt directory

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
		
