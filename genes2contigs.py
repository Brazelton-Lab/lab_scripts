#! /usr/bin/env python
# given list of gene ids (e.g. KEGG IDs), write table with each contig that has each gene id
# usage:
# python genes2contigs.py list-of-gene-ids.txt prokka.gff newfilename

# define files
import sys
gene_file = sys.argv[1]
gff_file = sys.argv[2]
newfilename = sys.argv[3]

# create list of gene ids
l = []
with open(gene_file) as genes:
	for gene in genes: 
		geneL = gene.split('\t')
		if len(geneL) > 1: gene = geneL[0]
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
			for contig in d[id]: newfile.write(contig + ',')	# for each contig that has this ID, write it to the file, separated by commas. last one has an orphan comma. sorry, deal with it
			newfile.write('\n')
		else: newfile.write(id + '\t' + '\n')		            
                        
