#! /usr/bin/env python
# outputs FASTA sequence corresponding to subject start and subject end in m8 file
# if subject has match over E value cutoff
# some subjects will be written multiple times, depending on how many hits they have

import sys
fastafilename = sys.argv[1]
blastfilename = sys.argv[2]

outfilename = blastfilename + '.besthits.fasta'
outfile = open(outfilename, 'a')

blast_file = open(blastfilename)
query_list = []
subject_list = []
count = 0
for line in blast_file:
	columns = line.split('\t')
	query = columns[0]
	subject = columns[1]
	Evalue = columns[10]
	try:
		Evalue = Evalue.split('e-')
		Evalue = Evalue[1]
		if int(Evalue) > 4:					# this is where you can change the evalue					
			count = count + 1
			query_list.append(query)
			subject_list.append(subject)
	except:
		if Evalue == '0':
			query_list.append(query)
			subject_list.append(subject)
			
print count,
print ' total hits above E value cutoff'
query_set = set(query_list)
subject_set = set(subject_list)
print len(query_set),
print ' unique queries with hits above E value cutoff'
print len(subject_set),
print ' unique subjects with hits above E value cutoff'

l=[] # keep track of subjects already found
blast_file = open(blastfilename)
for line in blast_file:
	columns = line.split('\t')
	subject = columns[1]
	if subject in l: pass # don't look for subjects already found
	else:
		start = int(columns[8])
		end = int(columns[9])
		Evalue = columns[10]
		
		if Evalue == '0': Evalue = 99
		elif 'e-' in Evalue:
			Evalue = Evalue.split('e-')
			Evalue = Evalue[1]
		else: Evalue = 1		
		
		fasta_file = open(fastafilename)
		status = 'stop'
		if int(Evalue) > 4:
			for line in fasta_file:			
				if line[0] == '>':
					status = 'stop'
					header = line.split(' ')
					header = header[0].replace('>','')
					header = header.replace('\n','')
					if subject == header:
						print subject					
						outfile.write('>')
						outfile.write(header)
						outfile.write('\n')
						l.append(header)
						status = 'go'
				elif status == 'go':		
					outfile.write(line)			
			
outfile.close()					
