#! /usr/bin/env python
# sort sequences in standard.fasta according to bar code
# files saved with filenames according to barcodes_distribution.txt

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
