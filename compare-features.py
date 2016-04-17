#! /usr/bin/env python
# given list of IDs (e.g. KEGG IDs), will write file of abundance of each ID in each function.tsv file in current directory
# usage:
# python compare-features.py list-of-IDs.txt newfilename

import sys
ids_file = sys.argv[1]
newfilename = sys.argv[2]

l = []
with open(ids_file) as ids:
	for id in ids: l.append(id.strip('\n'))

f = []
d = {}	
import glob
for filename in glob.glob('*.function.tsv'):
	with open(filename) as file:
		f.append(filename)
		for line in file:
			cols = line.split('\t')
			if cols[0] in l: 
				if cols[0] in d: d[cols[0]].append(cols[1].strip('\n'))
				else: d[cols[0]] = list(tuple([cols[1].strip('\n')]))	
			else: pass				
			
with open(newfilename,'w') as newfile:
	newfile.write('\t')
	for filename in f:
		filename = filename.rstrip('.function.tsv')
		newfile.write(filename + '\t')					# write header row
	newfile.write('\n')	
		
	for id in d:
		newfile.write(id + '\t')
		for value in d[id]: newfile.write(str(value) + '\t')	# assumes same order in d as in list f
		newfile.write('\n')
				
		
