# this script can be run before or after kraken or other programs
# only requirement is that contig name is in first column of input file and that the appropriate phylosift sequence_taxa_summary file is provided
# example usage:
# python table-phylosift.py file.length-GC-kraken.txt sequence_taxa_summary.txt

import sys
table_filename = sys.argv[1]
summary_filename = sys.argv[2]

outfilename = ''
words = table_filename.split('.')
for word in words[:-1]:
	outfilename = outfilename + word + '.' 
outfilename = outfilename + 'phylosift.txt'
print outfilename

startover = 'no'
with open(table_filename) as table_file:
	for line in table_file:
		genus = 'unknown'
		family = 'unknown'
		order = 'unknown'
		clas = 'unknown'
		
		cols = line.split('\t')
		name = cols[0]

		with open(summary_filename) as summary_file:
			for row in summary_file:
				fields = row.split('\t')
				contig = fields[0].split(' ')
				if name	== contig[0]:
					startover= 'yes'
					if fields[3] == 'genus':
						if float(fields[5]) > 0.6: genus = fields[4]
					if fields[3] == 'family':
						if float(fields[5]) > 0.6: family = fields[4]	
					if fields[3] == 'order':
						if float(fields[5]) > 0.6: order = fields[4]	
					if fields[3] == 'class':
						if float(fields[5]) > 0.6: clas = fields[4]	
				elif startover == 'no': pass
				else: 
					startover = 'no'
					break
								
		with open(outfilename,'a') as outfile: 
			outfile.write(line.replace('\n','\t'))
			outfile.write(genus + '\t' + family + '\t' + order + '\t' + clas + '\n')	
