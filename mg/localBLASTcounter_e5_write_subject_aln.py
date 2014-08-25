# outputs FASTA sequence corresponding to subject start and subject end in m8 file
# if subject has match over E value cutoff
# some subjects will be written multiple times, depending on how many hits they have

fasta_file = open('LCcontigs.fasta').read()
FASTAs = fasta_file.split('>')

print 'finished reading FASTA file'

outfile = open('CoxL_vs_LCcontigs_subjalns_e5.fasta', 'a')

blast_file = open('coxL_vs_LCcontigs.tblastn')
for line in blast_file.xreadlines():
	columns = line.split('\t')
	subject = columns[1]
	start = int(columns[8])
	end = int(columns[9])
	Evalue = columns[10]
	
	if Evalue == '0': Evalue = 99
	elif 'e-' in Evalue:
		Evalue = Evalue.split('e-')
		Evalue = Evalue[1]
	else: Evalue = 1		
	
	if int(Evalue) > 4:
		for FASTA in FASTAs:
			FASTA_split = FASTA.split('\n')
			header = FASTA_split[0]
			header = header.split(' ')
			header = header[0]
			#print subject,
			#print header
			if subject == header:
				#print subject
				seq = ''
				for remaining in FASTA_split[1:]:
					seq = seq + remaining
				first = start - 1
				last = end
				if last < first:
					outfile.write('>')
					outfile.write('c_')
					outfile.write(header)
					outfile.write('\n')
					for base in seq[first:last:-1]:
						if base == '\n': pass
						elif base == ' ': pass
						elif base == '\t': pass
						elif base.upper() == 'A': base = 'T'
						elif base.upper() == 'C': base = 'G'
						elif base.upper() == 'G': base = 'C'
						elif base.upper() == 'T': base = 'A'
						else: base = 'N'
						outfile.write(base)	
				elif first < last:
					outfile.write('>')
					outfile.write(header)
					outfile.write('\n')
					outfile.write(seq[first:last])
				else: print 'error'						
				outfile.write('\n')
				break
outfile.close()					