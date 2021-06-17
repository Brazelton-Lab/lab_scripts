import glob
import sys
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='writes a list of FASTA headers')
parser.add_argument('-e','--extension', required=True, help='file extension of each FASTA file (e.g. fa or fasta)')
parser.add_argument('-t','--table', required=False, help='write one table of results instead of a separate file for each input file (optional)')
args = parser.parse_args()

path = '*.' + args.extension.strip('.') # in case the user already provided the dot

if args.table:
	with open(args.table, 'w') as t:
		for fastafile in glob.glob(path):
			with open(fastafile) as f:
				for each in SeqIO.parse(f, "fasta"):
					t.write(fastafile)
					t.write('\t')
					t.write(each.id)
					t.write('\n')

else:
	for fastafile in glob.glob(path):
		outfilename = fastafile.replace(args.extension.strip('.'),'headers.txt')
		with open(fastafile) as f, open(outfilename, 'w') as o:
			for each in SeqIO.parse(f, "fasta"):
				o.write(str(each.id))
				o.write('\n')