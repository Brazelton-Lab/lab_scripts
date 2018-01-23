#! /usr/bin/env python

"""
sort sequences in standard.fasta according to bar code
files saved with filenames according to barcodes_distribution.txt

Copyright:

    make-group-file-with-barcodes  sort sequences in standard.fasta according to bar code

    Copyright (C) 2016  William Brazelton <comma-separated list of authors>

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
barcodes = open('barcodes_distribution_JGI1-2.txt')
for line in barcodes:
	if line[0] == "=": pass
	else: 
		columns = line.split('\t')
		sample = columns[0]
		barcode = columns[1]
		D[barcode] = sample
barcodes.close()

outfilename = infilename + '.group'
outfile = open(outfilename,'a')
		
fastas = open(infilename)
for line in fastas:
	#print line
	if line[0] == '>':
		header = line.strip('\n')
		words = line.split('#')
		fasta_barcode = words[1].strip('\n')
		fasta_sample = D[fasta_barcode]
		outfile.write(header.replace('>',''))
		outfile.write('\t')
		outfile.write(fasta_sample)
		outfile.write('\n')
	else: pass
fastas.close()	
outfile.close()
