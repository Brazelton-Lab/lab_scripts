# count A,C,T,G in a file ignoring headers designated as >
# write name of file and GC proportion to outfile

import sys
infilename = sys.argv[1]
outfilename = sys.argv[2]

infile_handle = open(infilename)
infile = infile_handle.read()
outfile = open(outfilename, 'a')
outfile.write('contig')
outfile.write('\t')
outfile.write('GC')
outfile.write('\t')
outfile.write('length')
outfile.write('\n')
		

infile = infile.split('>')
print (len(infile)-1),
print 'FASTA entries'

for FASTA in infile[1:]:
	A=0
	a=0
	C=0
	c=0
	G=0
	g=0
	T=0
	t=0
	FASTA_split = FASTA.split('\n')
	header = FASTA_split[0]
	seq = ''
	for remaining in FASTA_split[1:]:
		seq = seq + remaining
	for base in seq:
		if base == 'A':
			A=A+1
		if base == 'C':
			C=C+1
		if base == 'T':
			T=T+1
		if base == 'G':
			G=G+1
		if base == 'a':
			a=a+1
		if base == 'c':
			c=c+1
		if base == 'g':
			g=g+1
		if base == 't':
			t=t+1	
	GC = G + g + C + c
	AT = A + a + T + t
	ratio = float(GC) / (float(AT) + float(GC))
	
	outfile.write(header)
	outfile.write('\t')
	outfile.write(str(ratio))
	outfile.write('\t')
	outfile.write(str(len(seq)))
	outfile.write('\n')
infile_handle.close()
outfile.close()
		