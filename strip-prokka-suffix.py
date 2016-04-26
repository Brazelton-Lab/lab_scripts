#! /usr/bin/env python
# strip suffixes added to gene names in .gff files by prokka
# usage:
# strip-prokka-suffix.py 
# will run on any .gff file in same directory

import glob
for filename in glob.glob('*.gff'):
	with open(filename) as gff:
		newfilename = filename.replace('.gff','.good.gff')
		print 'writing' + newfilename
		with open(newfilename,'w') as newfile:
			for line in gff:
				if 'Name=' in line:
					fields = line.split(';')
					for field in fields:
						if 'Name=' in field:
							if '_' in field:
								name = field.split('_')
								newfile.write(name[0]+';')
						elif 'gene=' in field:
							if '_' in field:
								name = field.split('_')
								newfile.write(name[0]+';')
						elif if '\n' in field: newfile.write(field)		
						else: newfile.write(field+';')		
				else: newfile.write(line)				
