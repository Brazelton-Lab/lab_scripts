# keep rows in count table that match the OTU names in provided list
# assumes comma-delimited format
# only looks at first column of list file and first column of count table
# unaffected by additional columns, e.g. taxonomy columns are ok
# usage:
# python filter_count_table_by_otu.py count_table_filename list_filename

import sys
tablefilename = sys.argv[1]
listfilename = sys.argv[2]
outfilename = tablefilename.replace('.csv','.keep.csv')

otus = []

with open(outfilename, 'w') as outfile:
	with open(listfilename) as listfile:
		for row in listfile:
			cols = row.split(',')
			otus.append(cols[0].strip('\n'))
	with open(tablefilename) as tablefile:
		for row in tablefile:	# write the first line
			outfile.write(row)
			break
		for row in tablefile:
			cols = row.split(',')
			if cols[0] in otus: outfile.write(row)
			