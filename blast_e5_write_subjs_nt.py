#! /usr/bin/env python
# outputs FASTA sequence corresponding to subject start and subject end in m8 file
# if subject has match over E value cutoff
# some subjects will be written multiple times, depending on how many hits they have

"""
Copyright:

    blast_e5_wite_subjs_nt.py Output FASTA seqs from BLAST alignment
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
fastafilename = sys.argv[1]
blastfilename = sys.argv[2]

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
		if int(Evalue) > 4:
			count = count + 1
			#print Evalue
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
print 'extracting the longest nucleotide sequence from each subject that is aligned to a query....'
print 'total hits scanned: '

d={}
count = 0
blast_file = open(blastfilename)
for line in blast_file:
	columns = line.split('\t')
	subject = columns[1]
	length = columns[3]
	start = int(columns[8])
	end = int(columns[9])
	Evalue = columns[10]
	if Evalue == '0': Evalue = 99
	elif 'e-' in Evalue:
		Evalue = Evalue.split('e-')
		Evalue = Evalue[1]
	else: Evalue = 1
	if int(Evalue) > 4:
		count = count + 1
		print count,
		sys.stdout.flush()
		if subject not in d:
			d[subject] = ''
			proceed = 'yes'
		else:
			if int(len(d[subject])) < int(length):
				d[subject] = ''
				proceed = 'yes'
			else: proceed = 'no'

		if proceed == 'yes':
			fasta_file = open(fastafilename)
			status = 'stop'
			seq = ''
			for line in fasta_file:
				if line[0] == '>':
					status = 'stop'
					if seq == '': pass
					else:								# only proceed if a matching sequence was found in previous line
						seq = seq.replace('\n','')
						first = start - 1
						last = end
						if last < first:
							for base in seq[first:last:-1]:
								if base == '\n': pass
								elif base == ' ': pass
								elif base == '\t': pass
								elif base.upper() == 'A': base = 'T'
								elif base.upper() == 'C': base = 'G'
								elif base.upper() == 'G': base = 'C'
								elif base.upper() == 'T': base = 'A'
								else: base = 'N'
								d[subject] = d[subject] + base
						elif first < last:
							d[subject] = seq[first:last]					# write matching sequence from previous line before checking next FASTA sequence
						else: print 'error'
						break

					seq = ''
					header = line.split(' ')
					header = header[0].replace('>','')
					header = header.replace('\n','')
					if subject == header: status = 'go'

				elif status == 'go':		# only store sequence if matching header was found previously
					seq = seq + line

blast_file.close()
fasta_file.close()
outfilename = blastfilename + '.besthits.fasta'
outfile = open(outfilename, 'a')
# iterate through dictionary
for i in d:
	print i
	outfile.write('>')
	outfile.write(i)
	outfile.write('\n')
	outfile.write(d[i])
	outfile.write('\n')
print 'finished'
outfile.close()
