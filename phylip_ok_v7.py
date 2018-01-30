#! /usr/bin/env python
"""
rewrites FASTA file with only 10 characters in header so that it can be converted to phylip format.
assumes that the desired 10 characters are the last 10 characters in the element bound by the 1st and 2nd vertical lines"|"

usage:pythonphylip_ok.pyfilename

Copyright:

    phylip_ok_v7  rewrites FASTA file with only 10 characters in header so that it can be converted to phylip format.

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
filename=sys.argv[1]
file=open(filename)

outfilename=filename+'.phylipok'
outfile=open(outfilename,'a')

tablefilename = outfilename + '.remember.txt'
tablefile = open(tablefilename,'a')

count=0
for line in file:
	if line[0]=='>':
		count=count+1
		if 'MISEQ' in line:
			words=line.split('_')
			id=words[5]+'_'+words[6]
			id=id[:5]	#first5
			id=id.replace('>','')
			outfile.write('>')
			outfile.write(str(count)+'_')
			outfile.write(id)
			outfile.write('\n')
		elif'|'in line:
			words=line.split('|')
			if words[0] == '>gi':
				id = words[1]
				id = id[-5:]
			else:
				id=words[0]
				id=id[-5:]#last5
			id=id.replace('>','')
			outfile.write('>')
			outfile.write(str(count)+'_')
			outfile.write(id)
			outfile.write('\n')	
		else: 
			id = line[:6]	
			id=id.replace('>','')
			outfile.write('>')
			outfile.write(str(count)+'_')
			outfile.write(id)
			outfile.write('\n')
		tablefile.write(line.replace('\n',''))
		tablefile.write('\t')
		tablefile.write(str(count)+'_')
		tablefile.write(id)
		tablefile.write('\n')
	else: outfile.write(line)	
outfile.close()
tablefile.close()		
