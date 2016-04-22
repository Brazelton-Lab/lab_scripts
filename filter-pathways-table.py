#!/usr/bin/env python
# given output from pathways2contigs.py, filter table according to coverage or pathway
# usage:
# python filter-pathways-table.py file.IDs.contigs.tsv file.IDs.contigs.filter.tsv 0.5 all
# change 0.5 to the desired coverage threshold (inclusive, pathway with coverage of exactly 0.5 will be written to output file)
# change all to a specific pathway - only one allowed, but it can be general or specific
# i.e. this is valid input: Transporters;ABC_transporters,_prokaryotic_type;ABC-2_type_and_other_transporters;Antibiotic_ABC_transporter_[mD:m00248]
# and this is valid input: Transporters

import sys
filename = sys.argv[1]
newfilename = sys.argv[2]
coverage = sys.argv[3]
pathway = sys.argv[4]

if "'" in pathway: pass
else: print 'WARNING: pathway search is not accurate unless search term is enclosed in quotes'

l = 0
yes = 0
no = 0
with open(newfilename,'w') as newfile:
	with open(filename) as file:
		for line in file:
			l = l + 1
			cols = line.split('\t')
			if float(cols[1]) < float(coverage): no = no + 1
			else:
				if pathway == 'all': 
					newfile.write(line)
					yes = yes + 1
				elif pathway in line: 
					newfile.write(line)
					yes = yes + 1
				else: no = no + 1

print str(l) + ' lines in table'
print str(yes) + ' lines written to ' + newfilename
print str(no) + ' lines did not match criteria' 			

