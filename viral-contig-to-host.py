#! /usr/bin/env python3

# link viral contig to spacer (by creating a dictionary)
# link microbial contig to repeat (by creating a second dictionary)
# then link virus to putative host using the common CRISPR locus defined by spacer and repeat IDs in the dictionaries

import argparse

parser = argparse.ArgumentParser(description='writes a table of viral contigs linked to putative host contigs')
parser.add_argument('-s','--spacers', required=True, help='file containing viral contigs linked to spacers')
parser.add_argument('-r','--repeats', required=True, help='file containing host contigs linked to repeats')
parser.add_argument('-b','--bins', required=False, help='(optional) table of bin names and member contigs; format is two tab-delimited columns with bin name in first column and contig name in second column; if provided, bin names will be written to output file instead of contig names')
parser.add_argument('-o','--output', required=True, help='name of output file')
args = parser.parse_args()

Ds = {}		# dictionary of spacers and viral contigs
with open(args.spacers) as s:
	for line in s:
		cols = line.split(',')
		vcontig = cols[0]
		id = cols[26]
		if id[0] == 'G': 
			id = id.split('SP')
			G = id[0]
			Ds[vcontig] = G
		else: pass
		
Dr = {}		# dictionary of repeats and microbial contigs
with open(args.repeats) as r:
	for line in r:
		cols = line.split(',')
		mcontig = cols[0]
		id = cols[26]
		if id[0] == 'G': 
			id = id.split('DR')
			G = id[0]
			Dr[mcontig] = G
		else: pass

with open(args.output, 'w') as o:
	for vcontig in Ds:
		L = []
		for mcontig in Dr:
			if Dr[mcontig] == Ds[vcontig]: L.append(mcontig)
			else: pass
		o.write(vcontig)
		o.write(',')
		o.write(Ds[vcontig])
		for i in L:
			o.write(',')
			o.write(i)
		o.write('\n')
		
if args.bins:
	# create dictionary of bin names
	Db = {}		
	with open(args.bins) as b:
		for line in b:
			cols = line.split('\t')
			if cols[0] in Db:
				Db[cols[0]].append(cols[1].strip('\n'))
			else: 
				Db[cols[0]] = []
				Db[cols[0]].append(cols[1].strip('\n'))

	# write output file
	outfilename = args.output + '.bins.csv'
	with open(outfilename, 'w') as o:
		for vcontig in Ds:
			L = []
			for mcontig in Dr:
				if Dr[mcontig] == Ds[vcontig]: L.append(mcontig)
				else: pass
			
			B = []
			for i in L:
				for key, value in Db.items():
					if i in value: 
						B.append(key)
					
			o.write(vcontig)
			o.write(',')
			o.write(Ds[vcontig])
			for i in B:
				o.write(',')
				o.write(i)
			o.write('\n')
