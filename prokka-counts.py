#!
# reports summary counts from prokka gff file
# usage:
# python prokka-counts.py <sample.gff>

import sys
gff_file = sys.argv[1]

lines = 0
kegg = 0
uniprot = 0
pfam = 0
hamap = 0
clusters = 0
rrna = 0
trna = 0
other_rna = 0
other = 0
none = 0

with open(gff_file) as gff:
	for line in gff:
		if '##FASTA'in line: break
		if line[0] == '#': pass
		else: 
			lines = lines + 1
			if 'database' in line:
				fields = line.split(';')
				for field in fields:
					if field[:8] == 'database': 
						if 'kegg' in field: kegg = kegg + 1
						elif 'UniProtKB' in field: uniprot = uniprot + 1
						elif 'HAMAP' in field: hamap = hamap + 1
						elif 'CLUSTERS' in field: clusters = clusters + 1
						elif 'Pfam' in field: pfam = pfam + 1
						else: 
							print 'other: ' + str(field)
							other = other + 1
			elif 'barrnap' in line: rrna = rrna + 1
			elif 'Aragorn' in line: trna = trna + 1
			elif 'misc_RNA' in line: other_rna = other_rna + 1		
			else: none = none + 1

print str(lines) + '\t total lines in .gff'					
print str(kegg + uniprot + pfam + hamap + clusters + other + none) + '\t total protein-encoding genes' 
print str(rrna) + '\t rRNA'
print str(trna) + '\t tRNA'
print str(other_rna) + '\t other RNA'
print str(kegg) + '\t KEGG'
print str(uniprot) + '\t UniprotKB'
print str(pfam) + '\t Pfam' 
print str(hamap) + '\t HAMAP'
print str(clusters) + '\t CLUSTERS'
print str(other) + '\t other protein database'
print str(none) + '\t no annotation'
					
