# remove rows from count table with fewer than X total seqs
# assumes tab-delimited format
# assumes there IS a "total" column, as is standard with mothur count tables
# unaffected by additional columns, e.g. taxonomy columns are ok
# usage:
# python filter_count_table_mothur.py count_table_filename

# set threshold:
MINCOUNT = 2	# a MINCOUNT of 2 means total of 1. MINCOUNT of 4 means total of 2. because 'total' column causes this script's count to be twice the actual total


import sys
infilename = sys.argv[1]
if infilename[-3:] == 'tsv':
	outfilename = infilename.replace('.tsv','.nosingles.tsv')
else: outfilename = infilename + '.nosingles.tsv'

goodcount = 0
badcount = 0
with open(outfilename, 'w') as outfile:
	with open(infilename) as infile:
		for row in infile:	# write the first line
			outfile.write(row)
			break
		for row in infile:
			total = 0
			cols = row.split('\t')
			for i in cols:
				try: total = total + int(i)	#sums all numbers but ignores non-numbers
				except: pass
			if total > MINCOUNT: 
				goodcount = goodcount + 1
				outfile.write(row)
			else: badcount = badcount + 1

print 'Out of',
print (int(goodcount) + int(badcount)),
print 'total lines,',
print badcount,
print 'lines omitted in',
print outfilename