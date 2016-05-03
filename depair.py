#!/usr/bin/env python

'''seperate interleaves paired-end reads into two files

DEinterleave PAIRed-end reads (DEPAIR) accepts an interleaved FASTA or FASTQ
file containing interleaved paired end reads and outputs two FASTA or FASTQ
files, one containing the forward reads (R1) and the other containing the
reverse reads (R2).

Usage: depair.py [options]

    -fasta FASTA input file
    -fastq FASTQ input file
    -out output directory
    --version, -v prints version and exits
'''

__version__ = '0.10'

import argparse
import datetime
import os
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from bioinformatic_tools import qualityCheck

def deinterleave(in_file, out_file, file_type):
    #Sorts every in_file by alteranting sequences into two out_files
    with open(in_file, 'r') as in_handle:
        forward_file_name = out_file + '_R1.' + file_type
        with open(forward_file_name, 'a') as out_handle_forward:
            reverse_file_name = out_file + '_R2.' + file_type
            with open(reverse_file_name, 'a') as out_handle_reverse:
                counter = 0
                for seq_record in SeqIO.parse(in_handle, file_type):
                    if counter % 2 == 0:
                        SeqIO.write(seq_record, out_handle_forward,\
                                    file_type)
                    elif counter % 2 == 1:
                        SeqIO.write(seq_record, out_handle_reverse,\
                                    file_type)
                    counter += 1
    qualityCheck(out_file, file_type, in_file, forward_file_name,\
                  reverse_file_name)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Seperates a FASTA or ' + \
                                     'FASTQ file containing'+ \
                                     'interleaved paired-end reads into' + \
                                     'two seperate file')
    group = parser.add_mutually_exclusive_group(required = False)
    group.add_argument('-fasta',\
                       default = None,\
                        help = 'FASTA file to seperate')
    group.add_argument('-fastq',\
                       default = None,\
                       help = 'FASTQ file to seperate')
    parser.add_argument('-out',
                        default = None,\
                        help = 'name of output file, R1 and R2 added to name')
    parser.add_argument('--version', '-v',\
                        help = 'prints version and exit',
                        action = 'store_true')
    args = parser.parse_args()


    if args.version:
        print(__version__)
        sys.exit(0)
    elif args.out == None and args.fasta == None and args.fastq == None:
        print(__doc__)
        sys.exit(0)
    elif args.out == None:
        print('An output file name must be specified.')
        sys.exit(1)
    elif args.fasta != None:
        deinterleave(args.fasta, args.out, 'fasta')
    elif args.fastq != None:
        deinterleave(args.fastq, args.out, 'fastq')
    else:
        print('You must specify an input file.')
        sys.exit(1)
