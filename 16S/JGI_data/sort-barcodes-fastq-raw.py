#! /usr/bin/env python
# sort sequences in fastq according to bar code
# files saved with filenames according to barcodes_distribution.txt
# also sorts forward and reverse reads, so there will be two files per sample

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