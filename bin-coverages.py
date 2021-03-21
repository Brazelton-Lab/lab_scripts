#! /usr/bin/env python3

# sums (or calculates weighted average of) coverages for all contigs that belong to each bin provided by user
# bins defined by FASTA file containing contig names in FASTA header
# coverage table should be formatted like that produced by the table-coverage-from-bedtools.py script.

import sys
import os
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Sums coverages for all contigs belonging to each bin. If a lengths file is provided, then a weighted average will be calculated instead of a sum')

parser.add_argument('-i','--inputdir', required=True, help='directory containing FASTA files representing each bin')
parser.add_argument('-e','--extension', required=True, help='file extension of each FASTA file (e.g. fa or fasta)')
parser.add_argument('-c','--covfile', required=True, help='table of coverages where each row is one contig with the contig names in the first column')
parser.add_argument('-l', '--lengths', required=False, help='optional: tab-delimited table of lengths of each contig with the contig name in the first column and the length in the second column')
parser.add_argument('-o','--outfile', required=True, help='output table')
args = parser.parse_args()

path = args.inputdir + '/*.' + args.extension.strip('.') # in case the user already provided the dot

Dfinal = {} # final dictionary of bin coverages to be built during the operation of the script
Lfinal = {} # final dictionary of total contig lengths for each bin

# remember first line of cov table for later
with open(args.covfile) as c:
	for line in c:
		firstline = line
		break

# make dictionary of contig lengths if provided
if args.lengths:
	L = {}
	with open(args.lengths) as lf:
		for line in lf:
			cols = line.split('\t')
			name = cols[0]
			length = cols[1]
			L[name] = length
else: pass

import glob
for fastafile in glob.glob(path):
	
	# working with one bin at a time
	with open(fastafile) as f:
		contigList = []
		for each in SeqIO.parse(f, "fasta"):
			contigList.append(str(each.id))

	# construct dictionary of contig coverages for this bin
	d = {}
	# sum the lengths of all contigs in this bin
	total_length = 0

	with open(args.covfile) as c:
		for line in c:	break # skip first line
		for line in c:
			cols = line.split(',')
			id = cols[0]
			if id in contigList:
				newlist = []
				for i in cols[1:]:				# clean up the \n from the last entry and make all values into floats
					if isinstance(i,str): 
						if args.lengths:
							j = float(i.strip('\n')) * float(L[id])	# weight the average coverage by the length of the contig
							newlist.append(j)
						else: 
							j = float(i.strip('\n'))
							newlist.append(j)
					else: 
						if args.lengths:
							j = float(i) * float(L[id])
							newlist.append(float(j))
						else:
							j = float(i)
							newlist.append(float(j))
				d[id] = newlist
				if args.lengths: total_length = total_length + float(L[id])
	
	# sum coverages for all contigs in this bin
	l = []
	for i in d: l.append(d[i])	# make list of lists containing coverages
	m = [sum(x) for x in zip(*l)] 	# sum each column of the list
	
	bin_name = os.path.basename(fastafile)
	bin_name = bin_name.strip(args.extension)
	bin_name = bin_name.rstrip('.')
	Dfinal[bin_name] = m	# bin name, defined by fasta file, assigned to the new list of summed values
		
        # optional: store total length in a dictionary
	if args.lengths: Lfinal[bin_name] = total_length

# after Dfinal is completely built, write it to a file
with open(args.outfile, 'w') as o:
	try: firstline = firstline.replace('Contig','Bin')
	except: pass
	o.write(firstline)
	for i in Dfinal:
		o.write(i)
		for each in Dfinal[i]:
			o.write(',')
			if args.lengths: newcov = float(each) / float(Lfinal[i])
			else: newcov = each
			o.write(str(newcov))
		o.write('\n')
    		
    			
