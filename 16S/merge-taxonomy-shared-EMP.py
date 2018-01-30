#! /usr/bin/env python

"""
transfers taxonomy from .taxonomy file to each OTU in .shared.transposed file

to make shared.transposed file:
1)open .shared file in excel
2)copy and paste>transpose all data 
3)save as tab-delimited text
4)open in TextWrangler and change line breaks to Unix (LF)
5)then run script

usage:
python merge-taxonomy-shared.py myfile.taxonomy myfile.shared.transposed


Copyright:

    merge-taxonomy-shared-EMP  transfers taxonomy from .taxonomy file to each OTU in .shared.transposed file

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

import sys
taxfilename = sys.argv[1]
sharedfilename = sys.argv[2]
outfilename = sharedfilename + ".taxonomy"

D = {}
with open(taxfilename) as taxfile:
	for line in taxfile:
		cols = line.split('\t')
		otu = cols[0].split('_')
		#print cols
		otu = otu[-1]
		try: D[otu] = cols[1].replace(';','\t') 	
		except: pass
#print D
with open(sharedfilename) as sharedfile:
	for line in sharedfile:
		if "label" in line:
			with open(outfilename,'a') as outfile: outfile.write(line)
		elif "Group\t" in line: 
			with open(outfilename,'a') as outfile: outfile.write(line)
		elif "numOtus\t" in line: 
			with open(outfilename,'a') as outfile: outfile.write(line)
		else:
			with open(outfilename,'a') as outfile:
				outfile.write(line.strip('\n'))
				otu = line.split('\t')
				try: 
					otunum = otu[0]
					#print otunum
					outfile.write('\t' + D[otunum])			
				except: outfile.write(line)
