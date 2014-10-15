#! /usr/bin/env python
# edits mothur taxonomy summary file
# transfers last name that is not "unclassified" or "uncultured" to "unclassified" or "uncultured" assignment
# make sure that the file has default sorting (by rankID)

import sys
infilename = sys.argv[1]
outfilename = infilename + '.renamed.txt'
outfile = open(outfilename,'a')

infile = open(infilename)
for line in infile:
	if "unclassified" in line:
		columns = line.split('\t')
		tax = columns[2]
		newtax = tax + ' ' + lasttax
		outfile.write(columns[0])
		outfile.write('\t')
		outfile.write(columns[1])
		outfile.write('\t')
		outfile.write(newtax)
		for tab in columns[3:]:
			outfile.write('\t')
			outfile.write(tab)
	elif "uncultured" in line:
		columns = line.split('\t')
		tax = columns[2]
		newtax = tax + ' ' + lasttax
		outfile.write(columns[0])
		outfile.write('\t')
		outfile.write(columns[1])
		outfile.write('\t')
		outfile.write(newtax)
		for tab in columns[3:]:
			outfile.write('\t')
			outfile.write(tab)
	else: 
		outfile.write(line)
		columns = line.split('\t')
		lasttax = columns[2]
infile.close()
outfile.close()		
