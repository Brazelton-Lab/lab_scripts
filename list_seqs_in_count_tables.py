#! /usr/bin/env python

# makes list of OTUs found in one OR a second table
# ignores duplicates (i.e. only prints it once)
# assumes comma-delimited format
# usage:
# python list_seqs_in_count_tables.py count_table1 count_table2

import sys
table1 = sys.argv[1]
table2 = sys.argv[2]
outfilename = table1.replace('.csv','.list.csv')

otus = []

with open(table1) as t1:
	for row in t1:
		cols = row.split(',')
		otu = cols[0].strip('\n')
		otus.append(otu)
with open(table2) as t2:
	for row in t2:
		cols = row.split(',')
		otu = cols[0].strip('\n')
		if otu in otus: pass
		else: otus.append(otu)

with open(outfilename, 'w') as outfile:
	for otu in otus: outfile.write(otu +'\n')
			
