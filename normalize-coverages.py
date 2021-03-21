#! /usr/bin/env python3

# normalize coverage values according to read (or fragment) totals provided by user
# output is probably in units of TPM, assuming starting coverages are in reads (or fragments) per base or per kilobase and the totals are in the same units
# input is a table of coverages already normalized to length of target feature (e.g. the contig) produced, for example, by genomecov:
# bedtools genomecov -pc -ibam "$bam" > "$bam".cov

import sys
import argparse

parser = argparse.ArgumentParser(description='normalize coverage values according to read (or fragment) totals provided by user')
parser.add_argument('-c', '--cov', required=True, help='table of per-base coverages with sample names in first row and sequence names in first column')
parser.add_argument('-t', '--total', required=True, help='table of read or fragment totals for each sample listed in coverage table, with sample names in first column and read totals in second column')
parser.add_argument('-o', '--output', required=True, help='output filename')
parser.add_argument('-l', '--length_cutoff', help='optional: length cutoff such that only contigs of equal or greater length will be kept in the output')
parser.add_argument('-s', '--length_table', help='optional: lengths of contigs to be used for length_cutoff parameter, with contig names in first column and lengths in second column')
args = parser.parse_args()

# make dictionary of read totals for each sample
D = {}
with open(args.total) as t:
	for line in t:
		cols = line.split('\t')
		if len(cols) == 1: cols = line.split(',')       # tries to split by commas if no columns detected with tab delimiter
		D[cols[0]] = float(cols[1])	

if args.length_cutoff:		# if optional length_cutoff is provided
	# make list of acceptable contig names according to the length_cutoff
	C = []
	with open(args.length_table) as s:
		for line in s:
			cols = line.split('\t')
			if len(cols) == 1: cols = line.split(',')       # tries to split by commas if no columns detected with tab delimiter
			if int(cols[1]) >= int(args.length_cutoff): C.append(cols[0])
			else: pass

# index position of each sample name in cov table
I = []
with open(args.cov) as c:
	for line in c:
		cols = line.split('\t')
		if len(cols) == 1: cols = line.split(',')       # tries to split by commas if no columns detected with tab delimiter
		for i in cols: I.append((i).strip('\n'))
		break											# only the first line

# if optional length_cutoff is provided
if args.length_cutoff:	
	count = 0	# count number of contigs excluded from output file		
	with open(args.cov) as c, open(args.output, 'w') as o:
		for line in c:
			o.write(line)	# just the first line
			break
		for line in c:	# starting with second line
			cols = line.split('\t')
			if len(cols) == 1: cols = line.split(',')       # tries to split by commas if no columns detected with tab delimiter
			if cols[0] in C:								# checks if contig name is acceptable based on length_cutoff
				o.write(cols[0])
				for i in cols[1:]:
					cov = float(i)
					sample_pos = cols.index(i)
					sample = I[sample_pos]
					tot = D[sample]
					tpm = cov * (1 / tot) * 1000000				# equation borrowed from count.py in seq-annot package: https://github.com/Brazelton-Lab/seq-annot/blob/master/seq_annot/count.py
					o.write(',')
					o.write(str(tpm))
				o.write('\n')
			else: count = count + 1
	print(count),
	print("contigs excluded from output file")
	
# if optional length_cutoff is NOT provided
else:
	with open(args.cov) as c, open(args.output, 'w') as o:
		for line in c:
			o.write(line) 	# just the first line
			break
		for line in c:  # starting with second line
			cols = line.split('\t')
			if len(cols) == 1: cols = line.split(',')       # tries to split by commas if no columns detected with tab delimiter
			o.write(cols[0])
			for i in cols[1:]:
				cov = float(i)
				sample_pos = cols.index(i)
				sample = I[sample_pos]
				tot = D[sample]
				tpm = cov * (1 / tot) * 1000000				# equation borrowed from count.py in seq-annot package: https://github.com/Brazelton-Lab/seq-annot/blob/master/seq_annot/count.py
				o.write(',')
				o.write(str(tpm))
			o.write('\n')
