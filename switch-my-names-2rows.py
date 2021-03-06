#! /usr/bin/env python3

# this version expects the names to be organized by column in two long rows; there is a different version that expects the names to be organized by row in two long columns
# switches the names in the first column of provided file to those provided in a look-up "names" file
# names in the first ROW of the names file should match the names in the first column of the input file
# new file will contain the names provided in the second ROW of the names file
# output will be tab-delimited

import sys
import argparse

parser = argparse.ArgumentParser(description='switches your names according to your provided "names" file')
parser.add_argument('-i', '--input', help='input file, can be of any format as long as the first column contains the names to be switched')
parser.add_argument('-n', '--names', help='names file, should contain the original names in the first ROW and the new, desired names in the second ROW')
parser.add_argument('-o', '--output', help='name of the output file')
args = parser.parse_args()

l1 = []
l2 = []
D = {}
with open(args.names) as n:
	for line in n:
		cols = line.split('\t')			
		if len(cols) == 0: cols = line.split(',')	# tries to split by commas if no columns detected with tab delimiter
		for i in cols:
			l1.append(i.strip('\n'))
		break # only do it once
	for line in n: # do it for the second line
		cols = line.split('\t')			
		if len(cols) == 0: cols = line.split(',')	# tries to split by commas if no columns detected with tab delimiter
		for i in cols:
			if ',' in i:	# if there are multiple names separated by commas, take the first one
				j = i.split(',')
				i = j[0]
			l2.append(i.strip('\n'))
	D = dict(zip(l1, l2))	

with open(args.input) as i, open(args.output, 'w') as o:
	for line in i:
		cols = line.split('\t')
		if len(cols) == 0: cols = line.split(',')	# tries to split by commas if no columns detected with tab delimiter
		original = cols[0]
		
		try: 
			new = D[original]
			o.write(new)
			for i in cols[1:]:
				o.write('\t')	# change this if you want a different delimiter in output file
				o.write(i)

		except: print(original,"not found in names file")
		
		