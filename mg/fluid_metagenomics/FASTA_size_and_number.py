#! /usr/bin/env python

"""
report length of each FASTA sequence and report average

Copyright:

    FASTA_size_and_number  report length of each FASTA sequence and report average

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

infile = open('WHC2B_dereplicated_large_Geneious_contigs.fasta').read()
#outfile = open('length_FASTA_seq.txt', 'a')

length_list = []
fastas = infile.split('>')
for fasta in fastas[1:]:
	fasta_split = fasta.split('\n')
	header = fasta_split[0]
	seq = ""
	for line in fasta_split[1:]:	
		seq = seq + line	
	length = len(seq)
	length_list.append(length)
	#outfile.write(str(header) + '\t' + str(length) + '\n')

total = 0
for entry in length_list:
	total = int(total) + int(entry)
avg = float(total) / float(len(length_list))	
print len(length_list),
print " FASTA sequences"
print avg,
print " average length"
print min(length_list),
print " is shortest sequence"
print max(length_list),
print " is longest sequence"


	
