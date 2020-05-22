#! /usr/bin/env python
# export FASTA file of CDS protein sequences from a multi-record genbank file
# borrowed some from Cedar McKay's script (Rocap lab)

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import glob
for genbankfilename in glob.glob("*.gb*"):
	with open(genbankfilename) as g:
		if "LOCUS " in g.readline():
			fs = genbankfilename.split('.')
			fastafilename = fs[0]
			for e in fs[1:-1]:
				fastafilename = fastafilename + '.' + e
			fastafilename = fastafilename + '.faa'
			with open(fastafilename, 'w') as f:
				for seq_record in SeqIO.parse(g, "genbank"):
					definition =  seq_record.description
					i = 0
					for feature in seq_record.features:
						if feature.type=="CDS":
							i = i + 1
							id = definition + '_' + str(i)
							product = feature.qualifiers['product'][0]
			
							# stole location parsing from Cedar McKay
							#Go through some pain to make location human readable by adding 1 to first position
							location = str(int(str(feature.location.nofuzzy_start))+1) + ":" + str(feature.location.nofuzzy_end) 
							if feature.strand == 1:
								location = location + ' Forward'
							elif feature.strand == -1:
								location = location + ' Reverse'
							else:
								location = location + ' Could not determine strand'
			
							seq = Seq(feature.qualifiers['translation'][0])
							header = id + ' ' + product + ' ' + location
							newrecord = SeqRecord(seq, header, description='')
							SeqIO.write(newrecord,f,'fasta')
						else: pass
		else: pass
				
			
