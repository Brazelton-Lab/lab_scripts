# report length of each FASTA sequence and report average


import glob
import os
path = r'./'
for fileName in glob.glob(os.path.join(path, '*.fna')):
	infile = open(fileName)
	length_list = []
	length = 0
	for line in infile:
		if line[0] == '>': 
			if length == 0: pass
			else:
				length_list.append(length)
				length = 0
		else:
			length = length + len(line)
	length_list.append(length)		
				
	total = 0
	for entry in length_list:
		total = int(total) + int(entry)
	avg = float(total) / float(len(length_list))	
	print fileName
	print len(length_list),
	print " FASTA sequences"
	print total,
	print " total bps"
	print avg,
	print " average length"
	print min(length_list),
	print " is shortest sequence"
	print max(length_list),
	print " is longest sequence"


	