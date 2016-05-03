#! /usr/bin/env python
# count number of unique subjects with BLAST hits with greater than X E-value
# use m8 flag formatted output from BLAST
# results will print to screen

from sets import Set
import glob
import os
path = r'./'
for fileName in glob.glob(os.path.join(path, '*.tblastn')):
	print fileName
	infile = open(fileName)
	#infile_read = infile.read()
	#infile_rows = infile_read.split('\n')
	#print len(infile_rows)
	query_list = []
	subject_list = []
	count = 0
	for line in infile.xreadlines():
		count = count + 1
		columns = line.split('\t')
		query = columns[0]
		subject = columns[1]
		Evalue = columns[10]
		try:
			Evalue = Evalue.split('e-')
			Evalue = Evalue[1]
			if int(Evalue) > 4:
				#print Evalue
				query_list.append(query)
				subject_list.append(subject)
		except:
			if Evalue == '0':
				query_list.append(query)
				subject_list.append(subject)
	print count,
	print ' total hits'
	query_set = Set(query_list)
	subject_set = Set(subject_list)
	print len(query_set),
	print ' unique queries with hits above E value cutoff'
	print len(subject_set),
	print ' unique subjects with hits above E value cutoff'
	infile.close()
	
	# write unique subjects to file
	outfile_name = fileName.lstrip('./')
	outfile_name = outfile_name.replace('.tblastn', '_unique_subjectsE5.txt')
	outfile = open(outfile_name, 'a')
	for subject in subject_set:
		outfile.write(subject)
		outfile.write('\n')