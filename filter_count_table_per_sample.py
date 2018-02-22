#! /usr/bin/env python
# remove rows from count table with fewer than X total seqs as the maximum in any given sample
# assumes comma-delimited format
# assumes there is NO "total" column
# unaffected by additional columns, e.g. taxonomy columns are ok
# usage:
# python filter_count_table_per_sample.py count_table_filename 20
# change the '20' above to another number for the threshold of counts per sample

import sys
MINCOUNT = sys.argv[2]
infilename = sys.argv[1]
outfilename = infilename.replace('.tsv','.chopped.tsv')

with open(outfilename, 'w') as outfile:
	with open(infilename) as infile:
		for row in infile:	# write the first line
			if row =='#\n': pass
			else: 
				outfile.write(row)
				break
		for row in infile:
			max = 0
			cols = row.split('\t')
			for i in cols[1:]:			# skip first column
				if int(i) > max: max = int(i)
				else: pass
			if max > MINCOUNT: outfile.write(row)
