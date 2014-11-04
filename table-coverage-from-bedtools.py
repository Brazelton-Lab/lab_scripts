#! /usr/bin/env python

# calculates avg coverage per contig from bedtools genomecov output
# and adds the avg coverage number to the last column of a provided tab-delimited file
# only requirement for the tab-delimited file is that the first column contains contig names that match the contig names in the bedtools-generated file
# example usage:
# bedtools genomecov -ibam sample.sort.bam > sample.sort.bam.cov
# python table-coverage-from-bedtools.py contigs.txt sample.sort.bam.cov

import sys
table_filename = sys.argv[1]
bedtools_filename = sys.argv[2]

outfilename = ''
words = table_filename.split('.')
for word in words[:-1]:
	outfilename = outfilename + word + '.' 
outfilename = outfilename + 'cov.txt'
print outfilename

D = {}
with open(bedtools_filename) as bedtools:
	for line in bedtools:
		cols = line.split('\t')
		contig_name = cols[0]
		cov = cols[1]
		bases = cols[2]
		length = cols[3]
		total_coverage = float(cov) * float(bases) / float(length) 
		if contig_name in D: D[contig_name] = D[contig_name] + total_coverage
		else: D[contig_name] = total_coverage
		
with open(table_filename) as table_file:
	for row in table_file:
		columns = row.split('\t')
		with open(outfilename,'a') as outfile: 		
			outfile.write(row.replace('\n','\t'))
			try: outfile.write(str(D[columns[0]]) + '\n')
			except: 
				outfile.write('\n')
				print 'missing', columns[0]
