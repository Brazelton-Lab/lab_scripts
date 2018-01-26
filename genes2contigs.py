#! /usr/bin/env python

"""
given list of gene ids (e.g. KEGG IDs), write table with each contig that has each gene id

usage:
python genes2contigs.py list-of-gene-ids.txt prokka.gff newfilename
optionally, can also find corresponding ESOM data point IDs if provided a .names file
in this case, usage:
python genes2contigs.py list-of-gene-ids.txt prokka.gff newfilename esom.names

Copyright:

    genes2contigs  write table with each contig that has each gene id

    Copyright (C) 2016  William Brazelton <comma-separated list of authors>

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

# define files
import sys
gene_file = sys.argv[1]
gff_file = sys.argv[2]
newfilename = sys.argv[3]
status = 'no'
try:
	esom_file = sys.argv[4]
	status = 'yes'
except: pass

# if esom.names file is provided, make dictionary with contig name as key and data point(s) as values
if status == 'yes':
	e = {}
	with open(esom_file) as esom:
		for line in esom:
			cols = line.split('\t')
			if len(cols) > 1:
				contig = cols[1]
				contig = contig.split('_')
				contig = contig[0] + '_' + contig[1]
				if contig in e: e[contig].append(cols[0])   		# add this data point to an existing list of data points as the value (for contigs that are broken up into multiple data points)
				else: e[contig] = list(tuple([cols[0]]))             # create a new key with this contig and its first data point as the value

# create list of gene ids
l = []
with open(gene_file) as genes:
	for gene in genes: 
		geneL = gene.split('\t')
		if len(geneL) > 1: gene = geneL[0]
		gene = gene.strip()
		l.append(gene.strip('\n'))

# create dictionary from prokka gff file with each KEGG ID as a key and a list of contigs that have that KEGG ID as the value
d = {}
with open(gff_file) as gff:
	for line in gff:
		if '##FASTA'in line: break
		if line[0] == '#': pass
		elif 'gene_id=' in line:
			contig_id = line.split('\t')
			contig_id = contig_id[0].replace('-','_')
			gene_id = line.split('gene_id=')
			gene_id = gene_id[1].split(';')
			gene_id = gene_id[0]
			if gene_id in d: d[gene_id].append(contig_id)   # add this contig_id to an existing list of contig_ids as the value
			else: d[gene_id] = list(tuple([contig_id]))             # create a new key with this kegg_id and its host contig as the value

# iterate through list of gene ids, writing each one to file with corresponding contig
with open(newfilename, 'w') as newfile:
	for id in l:
		if id in d:
			newfile.write(id + '\t')
			for contig in d[id]: 
				newfile.write(contig + ',')	# for each contig that has this ID, write it to the file, separated by commas. last one has an orphan comma. sorry, deal with it
			if status == 'yes':
				newfile.write('\t')
				for contig in d[id]:			# iterate through contigs again to find esom data points if esom.names file provided 
					if contig in e:
						for data_point in e[contig]: newfile.write(data_point + ',')
					else: newfile.write('not in ESOM,')
			newfile.write('\n')
		else: newfile.write(id + '\t' + '\n')		            
                        
