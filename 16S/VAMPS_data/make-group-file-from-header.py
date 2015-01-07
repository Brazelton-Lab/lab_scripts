#! /usr/bin/env python
# sort sequences in standard.fasta according to bar code
# files saved with filenames according to barcodes_distribution.txt

import sys
infilename = sys.argv[1]

outfilename = infilename + '.group'
outfile = open(outfilename,'a')
		
fastas = open(infilename)
for line in fastas:
	#print line
	if line[0] == '>':
		header = line.strip('\n')
		header = header.split(' ')
		header = header[0]
		words = line.split('|')
		fasta_sample = words[2].split('_')
		if len(fasta_sample) < 3: sample = fasta_sample
		else: sample = fasta_sample[1:-1]

		outfile.write(header.replace('>',''))
		outfile.write('\t')
		
		outfile.write(sample[0])
		for i in sample[1:]:
			outfile.write('_')
			outfile.write(i)
		outfile.write('\n')
		
	else: pass
fastas.close()	
outfile.close()
