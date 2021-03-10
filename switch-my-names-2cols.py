#! /usr/bin/env python3

# switches the names in the first column of provided file to those provided in a look-up "names" file
# names in the first column of both provided files must match
# new file will contain the names provided in the second column of the provided names file
# output will be tab delimited

import sys
import argparse

parser = argparse.ArgumentParser(description='switches your names according to your provided "names" file')
parser.add_argument('-i', '--input', help='input file, can be of any format as long as the first column contains the names to be switched')
parser.add_argument('-n', '--names', help='names file, should contain the original names in the first column and the new, desired names in the second column')
parser.add_argument('-o', '--output', help='name of the output file')
args = parser.parse_args()

D = {}
with open(args.names) as n:
	for line in n:
		cols = line.split('\t')
		if len(cols) == 0: cols = line.split(',')	# tries to split by commas if no columns detected with tab delimiter
		D[cols[0]] = cols[1].strip('\n')	

with open(args.input) as i, open(args.output) as o:
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
		
		