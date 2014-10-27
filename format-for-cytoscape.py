# adds "X. to the front and ." to the end of each name in the first column if it begins with Otu
# used to format network attributes file for Cytoscape so that Otu names match
# assumes comma separated values and unix linebreaks
# usage:
# python format-for-cytoscape.py filename.csv

import sys
filename = sys.argv[1]
outfilename = filename.replace('.csv','.for-cys.csv')

status = 'no'
with open(filename) as file:
	for line in file:
		columns = line.split(',')
		if 'Otu' in columns[0]:
			status = 'yes'
			newname = '"X.' + columns[0] + '."'
			with open(outfilename,'a') as outfile:
				outfile.write(newname)
				for col in columns[1:]:
					outfile.write(',' + col)		
		else: 
			with open(outfilename,'a') as outfile: 
				outfile.write(line)	

if status == 'no': print 'Nothing happened. Perhaps your input file needs to be converted to Unix linebreaks? Delete the output file and try again.'				