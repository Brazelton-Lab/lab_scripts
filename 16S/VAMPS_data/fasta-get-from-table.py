#! /usr/bin/env python

"""

Copyright:

    fasta-get-from-table  export FASTA according to taxonomy listed in VAMPS spreadsheet

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

keep=[]

import csv
with open('DCO_BRZ_Bv4v5_TaxBySeq2.tsv','rb') as csvfile:
	reader = csv.reader(csvfile, delimiter='\t')
	for row in reader:
		if "Hydrogenophaga" in row:
			for column in row: 
				if len(column) > 300: keep.append(column)

#print keep
print len(keep),
print "Hydrogenophaga sequences identified in table"		

count = 0
outfilename = 'DCO_BRZ_Bv4v5_Hydrogenophaga2.fasta'
outfile = open(outfilename, 'a')
from Bio import SeqIO
for seq in SeqIO.parse("DCO_BRZ_Bv4v5.fa","fasta"):
	#print seq.seq
	if str(seq.seq) in keep: 
		count = count + 1
		SeqIO.write(seq,outfile,"fasta")	
		
print count,
print "Hydrogenophaga sequences exported to new FASTA file"		
		
outfile.close()	
