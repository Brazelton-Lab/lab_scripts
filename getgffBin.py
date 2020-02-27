#! /usr/bin/env python
# written by Julia McGonigle
# modified by Billy Brazelton to make it faster

import sys
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Add new contigs to faa bin')

parser.add_argument('-f','--fastabin', help='faa file to subset gff file from', required=True)
parser.add_argument('-g','--gffFile', help='annotated gff file to subset', required=True)
parser.add_argument('-o','--outFile', help='output gff file', required=True)
args = parser.parse_args()
    
contigList = []
with open(args.fastabin) as fastaFile:
    for each in SeqIO.parse(fastaFile, "fasta"):
        contigList.append(str(each.id))

contigList.sort()

with open(args.gffFile) as gff:
    with open(args.outFile, 'w') as outFile:
        for line in gff:
            cols = line.split('\t')
            id = cols[0]
	    if id in contigList: outFile.write(line)
            else: pass
