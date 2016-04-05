#!
# given a pathcoverage.tsv file from humann2 and a .gff file from Prokka, 
# write a table that lists pathway - coverage - KEGG ID - contig name which contains that KEGG ID
# also requires ontology file. for example: /srv/databases/ontology/foam
# usage:
# python pathways2table.py pathcoverage.tsv prokka.gff /srv/databases/ontology/foam

# define input files
import sys
pathway_file = sys.argv[1]
gff_file = sys.argv[2]
ontology_file = sys.argv[3]

# create name of new file
root = pathway_file.split('.')
root = root[0].split('/')
root = root[-1]
newfile = root + '.IDs.contigs.tsv'

# check that files are in proper order in command line
if pathway_file[-4:] == '.tsv': pass
if gff_file[-4:] == '.gff': pass
else:
	print 'Please enter files in this order: ESOM_bin.fasta assembly.gff assembly.faa assembly.ffn'
	sys.exit()

try: open(ontology_file,'r')
except: 
	print 'Unable to read ontology file. Did you provide it, or is it a permissions issue?'
	sys.exit()

# create dictionary with ontology file
o = {}
with open(ontology_file) as ontology:
	for line in ontology:
		line = line.strip('\n')
		cols = line.split('\t')
		o[cols[0].replace(' ','_')] = cols[1:]

# create dictionary from prokka gff file with each KEGG ID as a key and a list of contigs that have that KEGG ID as the value 
p = {}
with open(gff_file) as gff:
	for line in gff:
		if '##FASTA'in line: break
		if line[0] == '#': pass
		elif 'database=kegg' in line:
			contig_id = line.split('\t')
			contig_id = contig_id[0].replace('-','_')
			gene_id = line.split('gene_id=')
			gene_id = gene_id[1].split(';')
			kegg_id = gene_id[0]
			if kegg_id in p: p[kegg_id].append(contig_id)	# add this contig_id to an existing list of contig_ids as the value
			else: p[kegg_id] = list(tuple([contig_id]))		# create a new key with this kegg_id and its host contig as the value

# go through pathcoverage table, find component KEGG IDs for each pathway, then write the contig names that have those KEGG IDs
with open(newfile,'w') as new:
	with open(pathway_file) as pathways:
		for line in pathways:
			if 'Pathway' in line: pass
			elif 'UNMAPPED' in line: pass
			elif 'UNINTEGRATED' in line: pass
			else:
				cols = line.split('\t')
				for id in o[cols[0]]:	# cols[0] is the name of the pathway, matches ontology file
					if id in p:	# ignore KEGG IDs that were not in the Prokka gff file
						new.write(cols[0] + '\t' + cols[1].strip('\n') + '\t' + id + '\t') 
						for contig in p[id]: new.write(contig + ',')	# for each contig that has this KEGG ID, write it to the file, separated by commas. last one has an orphan comma. sorry, deal with it
        					new.write('\n')		
