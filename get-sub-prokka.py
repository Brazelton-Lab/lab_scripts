#!
# given a .fasta file of contigs, a .gff file from prokka, a .faa file from prokka, and a .ffn file from prokka,
# get entries in .faa and .ffn files that have Prokka IDs matching the contigs in the .fasta file as specified in the .gff file
# to be used to extract predicted genes and proteins from Prokka for a subset of the contigs fed to Prokka. For example, one ESOM bin from a larger assembly
# usage:
# python get-sub-prokka.py ESOM_bin.fasta assembly.gff assembly.faa assembly.ffn

import sys
contigfile = sys.argv[1]
gff_file = sys.argv[2]
faa_file = sys.argv[3]
ffn_file = sys.argv[4]

if gff_file[-4:] == '.gff': pass
if faa_file[-4:] == '.faa': pass
if ffn_file[-4:] == '.ffn': pass
else: 
	print 'Please enter files in this order: ESOM_bin.fasta assembly.gff assembly.faa assembly.ffn'
	sys.exit()

root = contigfile.split('.')
root = root[0].split('/')
root = root[-1]
newfaafile = root + '.faa'
newffnfile = root + '.ffn'

c = []
with open(contigfile) as contigs:
	for line in contigs: 
		if line[0] == '>': 
			header = line.strip('>')
			header = header.replace('-','_')
			header = header.replace(' ','_')
			header = header.split('_')
			header = header[0] + '_' + header[1]
			c.append(header.strip('\n'))

p = []
with open(gff_file) as gff:
	for line in gff:
		if '##FASTA'in line: break
		if line[0] == '#': pass
		else: 
			contig_id = line.split('\t')
			contig_id = contig_id[0].replace('-','_')
			if contig_id in c:
				prokka_id = line.split(';')
				if 'ID=' in line:
					prokka_id = prokka_id[0].split('ID=')
					p.append(prokka_id[1])

count = 0
from Bio import SeqIO
with open(newfaafile,'w') as newfaa:
	for faa in SeqIO.parse(faa_file,'fasta'):
		if faa.id in p: 
			SeqIO.write(faa,newfaa,'fasta')
			count = count + 1

print count,
print 'entries from .faa file written to',
print newfaafile
count = 0
	
from Bio import SeqIO
with open(newffnfile,'w') as newffn:
	for ffn in SeqIO.parse(ffn_file,'fasta'):
		if ffn.id in p: 
			SeqIO.write(ffn,newffn,'fasta')
			count = count + 1

print count,
print 'entries from .ffn file written to',
print newffnfile
