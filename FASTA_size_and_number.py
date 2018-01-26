#! /usr/bin/env python


"""
report length of each FASTA sequence and report average
usage: python FASTA_size_and_number.py infilename
usage: python FASTA_size_and_number.py infilename report sort
second usage will report name and length of each FASTA sequence
leave out "sort" if you don't want to change the order in the file

Copyright:

    FASTA_size_and_number.py  Report length of each sequence and average len
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
infilename = sys.argv[1]

infile = open(infilename).read()
#outfile = open('length_FASTA_seq.txt', 'a')

length_list = []
header_list = []
fastas = infile.split('>')
for fasta in fastas[1:]:
	fasta_split = fasta.split('\n')
	header = fasta_split[0]
	header_list.append(header)
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

zipped = zip(length_list, header_list)
s = sorted(zipped)

try:
	report = sys.argv[2]
	try:
		issort = sys.argv[3]
		for i in s:
			print i
	except:
		for i in zipped:
			print i
except: pass




