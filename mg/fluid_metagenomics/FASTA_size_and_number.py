# report length of each FASTA sequence and report average

infile = open('WHC2B_dereplicated_large_Geneious_contigs.fasta').read()
#outfile = open('length_FASTA_seq.txt', 'a')

length_list = []
fastas = infile.split('>')
for fasta in fastas[1:]:
	fasta_split = fasta.split('\n')
	header = fasta_split[0]
	seq = ""
	for line in fasta_split[1:]:	
		seq = seq + line	
	length = len(seq)
	length_list.append(length)
	#outfile.write(str(header) + '\t' + str(length) + '\n')

total = 0
for entry in length_list:
	total = int(total) + int(entry)
avg = float(total) / float(len(length_list))	
print len(length_list),
print " FASTA sequences"
print avg,
print " average length"
print min(length_list),
print " is shortest sequence"
print max(length_list),
print " is longest sequence"


	