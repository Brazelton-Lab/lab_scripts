# export FASTA according to taxonomy listed in VAMPS spreadsheet

keep=[]

import csv
with open('DCO_BRZ_Bv4v5_TaxBySeq2.tsv','rb') as csvfile:
	reader = csv.reader(csvfile, delimiter='\t')
	for row in reader:
		if "Hydrogenophaga" in row:
			for column in row: 
				if len(column) > 300: keep.append(column)

#print keep
print len(keep),
print "Hydrogenophaga sequences identified in table"		

count = 0
outfilename = 'DCO_BRZ_Bv4v5_Hydrogenophaga2.fasta'
outfile = open(outfilename, 'a')
from Bio import SeqIO
for seq in SeqIO.parse("DCO_BRZ_Bv4v5.fa","fasta"):
	#print seq.seq
	if str(seq.seq) in keep: 
		count = count + 1
		SeqIO.write(seq,outfile,"fasta")	
		
print count,
print "Hydrogenophaga sequences exported to new FASTA file"		
		
outfile.close()	