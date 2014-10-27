# parses p values output from rcor.test in R
# usage:
# python parse-pvalues.py yourfilename.rcortest.pvalues

import sys
filename = 'myfilename.rcortest.pvalues.csv'
f = open(filename)

outfilename = filename.replace('csv','txt')
outfile = open(outfilename, 'a')

for line in f:
	columns = line.split(',')
	if columns[0] == '""': pass	# skips first line
	elif len(columns) < 2: pass	# skips last line
	else:	
		pvalue = columns[3]
		if pvalue == '"NA\n"': print 'Your input data probably includes columns containing only zeroes. You should delete these and repeat the R commands.'
		else: outfile.write(pvalue)
 
f.close()
outfile.close() 
 