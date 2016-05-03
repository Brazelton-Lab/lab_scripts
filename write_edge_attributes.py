#! /usr/bin/env python

# THIS SCRIPT NOT NECESSARY - JUST CLICK ON THE COLUMN IN CYTOSCAPE IMPORT WINDOW THAT YOU WANT TO BE THE EDGE ATTRIBUTES

# from file with variable names and correlation coefficient, export Edge Attributes file for importing with Cytoscape
# infile must end with ".txt"
# usage: python write_edge_attributes.py filename

import sys
infilename = sys.argv[1]
outfilename = infilename.replace('.txt','_edge_attr.txt')
infile = open(infilename)
outfile = open(outfilename,'a')
outfile.write('Continuous edge Attr1\n')

for line in infile.xreadlines():
	if line[:4] == "var1": pass
	elif len(line) < 2: pass
	else:
		columns = line.split('\t')
		outfile.write(columns[0])
		if float(columns[2]) > 0: outfile.write(' (pos) ')
		elif float(columns[2]) < 0: outfile.write(' (neg) ')
		elif float(columns[2]) == 0: outfile.write(' (neutral) ')
		else: outfile.write(' (error) ')
		outfile.write(columns[1])
		outfile.write(' = ')
		outfile.write(columns[2])
		outfile.write('\n')
infile.close()
outfile.close()
