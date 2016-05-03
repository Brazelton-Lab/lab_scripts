#! /usr/bin/env python
# given list of IDs (e.g. KEGG IDs), will write file of abundance of each ID in each function.tsv file in current directory
# usage:
# python compare-features.py list-of-IDs.txt newfilename

import sys
ids_file = sys.argv[1]
newfilename = sys.argv[2]

l = []
with open(ids_file) as ids:
	for line in ids: 
		cols = line.split('\t')
		id = cols[0].strip()
		l.append(id.strip('\n'))

f = []	# need to make defined list of filenames to make sure header row in new file is in correct order
import glob
for filename in glob.glob('*.function.tsv'):
	f.append(filename)
f.sort()

d = {}	
for id in l: d[id] = []		# establish entry for each ID in the list

count = 0
for filename in f:
	count = count + 1
	with open(filename) as file:
		for line in file:
			cols = line.split('\t')
			if cols[0] in l: d[cols[0]].append(cols[1].strip('\n'))
			else: pass
	for id in d: 
		if len(d[id]) == count: pass
		else: d[id].append('0')
			
with open(newfilename,'w') as newfile:
	newfile.write('\t')
	for filename in f:
		filename = filename.rstrip('.function.tsv')
		newfile.write(filename + '\t')					# write header row
	newfile.write('\n')	
		
	for id in l:
		newfile.write(id + '\t')
		if id in d:
			for value in d[id]: newfile.write(str(value) + '\t')	# assumes same order in d as in list f
			newfile.write('\n')
		else: newfile.write('\n')
				
		
