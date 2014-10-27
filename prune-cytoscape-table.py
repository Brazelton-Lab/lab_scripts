# prune tab-delimited text file according to user-provided list of IDs
# i.e. will export those rows that include one of the user-provided IDs
# usage:
# python prune-cytoscape-table.py your_cors_and_pqvalues.05.txt filename-with-IDs.txt AND

import sys
filename1 = sys.argv[1]
filename2 = sys.argv[2]
try: logic = sys.argv[3]
except: 
	print "please type either 'AND' or 'OR' without the quotes"
	sys.exit()
outfilename = filename1.replace('txt',logic)
outfilename = outfilename + '.' + filename2

file2 = open(filename2)
IDs = []
for line in file2: 
	ID = line.strip('\n')
	IDs.append(ID)

outfile = open(outfilename, 'a')
file1 = open(filename1)
for line in file1:
	if line[:4] == 'var1': outfile.write(line)	# keep header row
	columns = line.split('\t')
	
	if logic == 'AND':
		if columns[0] in IDs: 
			if columns[1] in IDs: outfile.write(line)
	elif logic == 'IF':
		if columns[0] in IDs: outfile.write(line)
		elif columns[1] in IDs: outfile.write(line)
		else: pass
	else: print "please type either 'AND' or 'OR' without the quotes"
file1.close()
file2.close()
outfile.close()
