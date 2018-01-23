#! /usr/bin/env python

"""
add a column to an OTU table that contains information from the SourceTracker2 per-feature assignment files

usage:
srun sourcetracker_to_otu_table.py count_table.txt

requires OTU table of same format as that used for SourceTracker; i.e. header row begins with #OTU
will look for SourceTracker per-feature assignment files in same directory named *.feature_table.txt
sums the probabilities across all feature_table.txt files and reports the category with highest summed probability

will store data as a dictionary. each source category is a key in the dictionary, and its value is a dictionary assigning a list of probabilities to each OTU
so it looks like this:
'category1'['ACTGCAACGTGCAGTG'] = [3,400,2,5]
the list of probabilities is assembled across all feature_table files

Copyright:

    sourcetracker_to_otu_table  add a column to an OTU table

    Copyright (C) 2016  William Brazelton <comma-separated list of authors>

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

lf = []	# list of filenames
dd = {}	# dictionary of dictionaries, one for each source category in feature_tables

# build dd each feature_table file at a time
import glob
for filename in glob.glob('*.feature_table.txt'):
	lf.append(filename)
	with open(filename) as f:
		firstline = f.readline().strip()				# all OTUs are on the first line of the feature_table files
		otus = firstline.split('\t')					
		otus = ['blank'] + otus 						# first item in otus list should be blank so that first OTU is otus[1]
		for line in f:
			d = {}
			cols = line.split('\t')
			source = cols[0]
			if source in dd: pass						
			else: dd[source] = {}						# creates a dictionary for each source
			for i in cols[1:-1]:						# builds dictionary of dictionaries, each containing OTU name assigned to list of probabilities
				n = cols.index(i)
				if otus[n] in dd[source]: dd[source][otus[n]].append(int(i))		# extends list of probability values if this OTU has been seen before
				else: dd[source][otus[n]] = list(tuple([int(i)]))				# links name of OTU to probability value from same column in this line 

# obtain list of OTUs from OTU table
otus = []
import sys
count_table_filename = sys.argv[1]
with open(count_table_filename) as c:
	for row in c:
		cols = row.split('\t')
		if row[0] == '#': pass
		else: otus.append(cols[0])
		
# evaluate results: sum probability values for each OTU for each source across all feature_tables
dw = {} 		# dictionary of winning source for each OTU
for otu in otus:
	highest = 0
	winner = 'blank'
	for i in dd:
		if otu in dd[i]: 
			total = sum(dd[i][otu])
			if int(total) > int(highest): 
				winner = i
				highest = total	
	dw[otu] = winner

# define outfilename
words = count_table_filename.split('.')
outfilename = words[0]
for i in words[1:-1]:									
	outfilename = outfilename + '.' + i 
outfilename = outfilename + '.st.tsv'
print 'writing',
print outfilename

# open OTU table and write result from feature_table files in each line
with open(outfilename, 'w') as o:
	with open(count_table_filename) as c:
		for row in c:
			if row[:4] == '#OTU': 
				o.write(row.strip() + '\t' + 'SourceTracker' + '\n')	
			elif '#' in row: pass
			else:
				cols = row.split('\t')
				seq_name = cols[0]
				if seq_name in dw:
					o.write(row.strip())
					o.write('\t' + dw[seq_name] + '\n')
				else: o.write(row)
