#! /usr/bin/env python

"""
calculates avg coverage per contig from bedtools genomecov output
collects contig coverages from all .cov files in working directory into a single .csv file

Copyright:

    table-coverage-from-bedtools  calculates avg coverage per contig from bedtools genomecov output

    Copyright (C) 2020  William Brazelton

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
import glob
Dall = {}
print 'collecting names of all contigs......',
for filename in glob.glob("*.cov"):			# just collect all contig names by looping through all .cov files
	with open(filename) as bedtools:
		for line in bedtools:
			cols = line.split('\t')
			contig_name = cols[0]
			Dall[contig_name] = []	
# now Dall has complete list of contigs
print 'done'

print 'collecting contig coverages from samples:'
samples = []
for filename in glob.glob("*.cov"):
	samples.append(filename)
	print filename
	d = {}
	with open(filename) as bedtools:
		for line in bedtools:
			cols = line.split('\t')
			contig_name = cols[0]
			cov = cols[1]
			bases = cols[2]
			length = cols[3]
			total_coverage = float(cov) * float(bases) / float(length) 
			if contig_name in d: d[contig_name] = d[contig_name] + total_coverage
			else: d[contig_name] = total_coverage
		
		# loop through all contigs in Dall
		for contig in Dall:
			if contig in d: Dall[contig].append(d[contig]) # add the coverage for this sample to the list of coverages for this contig in Dall
			else: Dall[contig].append(0)					# if this contig was not in this sample, add a 0 to the list of coverages for this contig in Dall

print 'writing to contig-coverages.csv......',
outfilename = "contig-coverages.csv"
with open(outfilename,'w') as outfile:
	for sample in samples:
		outfile.write(',' + sample)
	outfile.write('\n')
	for contig in Dall:
		outfile.write(contig)
		for i in Dall[contig]:
			outfile.write(',' + str(i))
		outfile.write('\n')
print 'done'
