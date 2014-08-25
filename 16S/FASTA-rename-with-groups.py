# renames FASTA file according to .groups file
# use when sending FASTA file from mothur to Qiime
# usage:
# python FASTA-rename-with-groups.py filename.fa filename.group

import sys
fastafilename = sys.argv[1]
groupfilename = sys.argv[2]

outfilename = ''
splitname = fastafilename.split('.')
for word in splitname[:-1]:
	outfilename = outfilename + word + '.'
outfilename = outfilename + 'groups.fa'

D = {}
with open(groupfilename) as groupfile:
	for line in groupfile:
		columns = line.split('\t')
		D[columns[0]] = columns[1].strip('\n')

count = 0
with open(fastafilename) as fasta:
	for line in fasta:
		if '>' in line:
			count = count + 1
			name = line.strip('\n')
			name = name.replace('>','')
			try: 
				groupname = D[name]
				with open(outfilename,'a') as outfile: outfile.write('>' + groupname + '_' + str(count) + ' ' + name + '\n')
			except: print name + ' not found'
		else: 
			with open(outfilename,'a') as outfile: outfile.write(line)	

	
