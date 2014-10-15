#! /usr/bin/env python

# replaces names in tree file with original names as provided by .remember.txt table
# usage:
# remember-names.py file.tre file.remember.txt

import sys
treefilename = sys.argv[1]
rememberfilename = sys.argv[2]
outfilename = treefilename + '.newnames.tre'

d = {}
with open(rememberfilename) as rfile:
	for line in rfile:
		cols = line.split('\t')
		d[cols[1].strip('\n')] = cols[0].strip('>')	

import re
with open(treefilename) as treefile:
	for tree in treefile:
		names = tree.split(',')
		for name in names:
			if name in d: 
				name_no_pars = name.replace('(','')
				name_no_pars = name_no_pars.replace(')','')
				oldname = name_no_pars.replace('QUERY___','')
				newname = d[oldname]
				print oldname,
				print ' --> ',
				print newname
			elif ':' in name: 
				name_before_colon = name.split(':')
				name_no_pars = name_before_colon[0].replace('(','')
				name_no_pars = name_no_pars.replace(')','')
				oldname = name_no_pars.replace('QUERY___','')
				newname = d[oldname]
				print oldname,
				print ' --> ',
				print newname
			else: 
				print 'name not recognized'
				break
			# write everything in between the commas (i.e. 'name') but replace oldname with newname
			if '\n' in name: # i.e. the last one, don't want to write another comma, which will cause an error in reading the tree
				if 'QUERY___' in name:
					with open(outfilename,'a') as outfile: outfile.write(name.replace('QUERY___'+oldname,newname))	
				else:
					with open(outfilename,'a') as outfile: outfile.write(name.replace(oldname,newname))	
			else:
				if 'QUERY___' in name:
					with open(outfilename,'a') as outfile: outfile.write(name.replace('QUERY___'+oldname,newname) + ',')	
				else:
					with open(outfilename,'a') as outfile: outfile.write(name.replace(oldname,newname) + ',')	
