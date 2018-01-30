#! /usr/bin/env python

"""
sort sequences in fastq according to bar code
files saved with filenames according to barcodes_distribution.txt
also sorts forward and reverse reads, so there will be two files per sample

Copyright:

    sort-barcodes-fastq-raw  sort sequences in fastq according to bar code

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

D={}
barcodes = open('barcodes_distribution.txt')
for line in barcodes:
	if line[0] == "=": pass
	else: 
		columns = line.split('\t')
		sample = columns[0]
		barcode = columns[1]
		D[barcode] = sample
barcodes.close()
		
fastas = open(infilename)
for line in fastas:
	#print line
	if line[:6] == '@MISEQ':
		try: outfile.close()
		except: pass
		header = line.split(':')
		fasta_barcode = header[-1].strip('\n')
		try: fasta_sample = D[fasta_barcode]
		except: fasta_sample = 'notfound'
		
		direction = header[-4]
		direction = direction.split(' ')
		direction = direction[1]
		
		outfile = open(fasta_sample+'.'+direction+'.fastq','a')
		outfile.write(line)
	else: outfile.write(line)
fastas.close()	
