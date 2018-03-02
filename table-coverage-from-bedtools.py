#! /usr/bin/env python

"""
calculates avg coverage per contig from bedtools genomecov output
and adds the avg coverage number to the last column of a provided tab-delimited file
only requirement for the tab-delimited file is that the first column contains contig names that match the contig names in the bedtools-generated file
.fasta file can be provided instead of tab-delimited file if desired - this functionality is triggered by the filename ending in exactly '.fasta'

example usage:
bedtools genomecov -ibam sample.sort.bam > sample.sort.bam.cov
python table-coverage-from-bedtools.py contigs.txt sample.sort.bam.cov
or
python table-coverage-from-bedtools.py contigs.fasta sample.sort.bam.cov

Copyright:

    table-coverage-from-bedtools  calculates avg coverage per contig from bedtools genomecov output

    Copyright (C) 2016  William Brazelton

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import sys
table_filename = sys.argv[1]
bedtools_filename = sys.argv[2]

outfilename = ''
words = bedtools_filename.split('.')
for word in words[:-1]:
	outfilename = outfilename + word + '.' 
outfilename = outfilename + 'cov.txt'
print 'will write to', outfilename

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

status = 'none'
if table_filename[-3:] == '.fa': status = 'fasta'
elif table_filename[-4:] == '.fna': status = 'fasta'
elif table_filename[-6:] == '.fasta': status = 'fasta'
else: status = 'table'
if status == 'fasta':
	l = []
	from Bio import SeqIO
	for fasta in SeqIO.parse(table_filename,"fasta"):
		l.append(fasta.id)
		
	with open(outfilename,'w') as outfile:
		outfile.write('Contig' + '\t' + bedtools_filename + '\n')
		for contig in l:
			outfile.write(contig + '\t')
			try: outfile.write(str(D[contig]) + '\n')
			except:
				outfile.write(contig + '\t' + '\n')
				print 'missing', contig

elif status == 'table':		
	with open(outfilename,'w') as outfile:
		with open(table_filename) as table_file:
			first_line = table_file.readline()
			outfile.write(first_line.replace('\n','\t'))
			outfile.write(bedtools_filename + '\n')
			for row in table_file:
				columns = row.split('\t') 		
				outfile.write(row.replace('\n','\t'))
				try: outfile.write(str(D[columns[0]]) + '\n')
				except: 
					outfile.write('\n')
					print 'missing', columns[0]
