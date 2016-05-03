#! /usr/bin/env python

#rewritesFASTAfilewithonly10charactersinheadersothatitcanbeconvertedtophylipformat
#assumesthatthedesired10charactersarethelast10charactersintheelementboundbythe1stand2ndverticallines"|"
#usage:pythonphylip_ok.pyfilename

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
