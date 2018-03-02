#! /usr/bin/env python

"""
strip suffixes added to gene names in .gff files by prokka

usage:
strip-prokka-suffix.py 
will run on any .gff file in same directory

Copyright:

    strip-prokka-suffix  strip suffixes added to gene names in .gff files by prokka

    Copyright (C) 2016  William Brazelton

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

import glob
for filename in glob.glob('*.gff'):
	with open(filename) as gff:
		newfilename = filename.replace('.gff','.good.gff')
		print 'writing ' + newfilename
		with open(newfilename,'w') as newfile:
			for line in gff:
				if 'Name=' in line:
					fields = line.split(';')
					for field in fields:
						if 'Name=' in field:
							if '_' in field:
								name = field.split('_')
								newfile.write(name[0]+';')
							else: newfile.write(field+';')
						elif 'gene=' in field:
							if '_' in field:
								name = field.split('_')
								newfile.write(name[0]+';')
							else: newfile.write(field+';')
						elif '\n' in field: newfile.write(field)		
						else: newfile.write(field+';')		
				else: newfile.write(line)				
