#! /usr/bin/env python3

# sums coverages for all contigs that belong to each bin provided by user
# bins defined by FASTA file containing contig names in FASTA header
# coverage table should be formatted like that produced by the table-coverage-from-bedtools.py script.

import sys
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Sums coverages for all contigs belonging to each bin')

parser.add_argument('-i','--inputdir', help='directory containing FASTA files representing each bin', required=True)
parser.add_argument('-e','--extension', help='file extension of each FASTA file (e.g. fa or fasta)', required=True)
parser.add_argument('-c','--covfile', help='table of coverages where each row is one contig with the contig names in the first column', required=True)
parser.add_argument('-o','--outfile', help='output table', required=True)
args = parser.parse_args()

path = args.inputdir + '/*.' + args.extension.strip('.') # in case the user already provided the dot

Dfinal = {} # final dictionary of bin coverages to be built during the operation of the script

# remember first line of cov table for later
with open(args.covfile) as c:
	for line in c:
		firstline = line
		break

import glob
for fastafile in glob.glob(path):
	
	# working with one bin at a time
	with open(fastafile) as f:
		contigList = []
		for each in SeqIO.parse(f, "fasta"):
			contigList.append(str(each.id))

	# construct dictionary of contig coverages for this bin
	d = {}
	with open(args.covfile) as c:
		for line in c:	break # skip first line
		for line in c:
			cols = line.split(',')
			id = cols[0]
			if id in contigList:
				newlist = []
				for i in cols[1:]:				# clean up the \n from the last entry and make all values into floats
					if isinstance(i,str): 
						newlist.append(float(i.strip('\n')))
					else: newlist.append(float(i))
				d[id] = newlist
				
	# sum coverages for all contigs in this bin
	for contig in d:
		l = []
		for i in d: l.append(d[i])	# make list of lists containing coverages
		m = [sum(x) for x in zip(*l)] 	# sum each column of the list
		
		bin_name = fastafile.strip(args.extension)
		bin_name = bin_name.strip(args.inputdir)
		bin_name = bin_name.strip('/')
		Dfinal[bin_name.rstrip('.')] = m	# bin name, defined by fasta file, assigned to the new list of summed values

# after Dfinal is completely built, write it to a file
with open(args.outfile, 'w') as o:
	try: firstline = firstline.replace('Contig','Bin')
	except: pass
	o.write(firstline)
	for i in Dfinal:
		o.write(i)
		for each in Dfinal[i]:
			o.write(',')
			o.write(str(each))
		o.write('\n')
    		
    			