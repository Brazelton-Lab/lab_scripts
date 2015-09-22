#!/usr/bin/env python

'''uniqueFastaId alters the id of each FASTA entry in a file

Some programs such as IDBA-UD may create multiple fasta/q files where
the ID of the entries in each file are identical. uniqueFastaId simply takes a
fasta/q file and appends a user defined string to each Fasta Header ID to
ensure that each file will not have headers that clash if the files need to be
concatenated for use by various programs. Unqique adds a number to each ID to
avoid repeats.

Usage: uniqueFastaId.py <fasta/q file> <--string to append> <--unique> [options]

    --version, -v prints version and exits
'''

__version__ = '0.12'

import argparse
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from bioinformatic_tools import fastaOrFastq
from bioinformatic_tools import qualityCheck

def appendStringToHeaders(in_file, string_to_append):
    #Appends string to each FASTA Header ID
    fileType = fastaOrFastq(in_file)
    with open(in_file, 'r') as in_handle:
        out_file = ''
        for i in in_file.split('.')[0:-1]:
            out_file += i
        out_file += '.unique.' + fileType
        with open(out_file, 'a') as out_handle:
            for seq_record in SeqIO.parse(in_handle, fileType):
                seq_record.id =  '_' + string_to_append
                SeqIO.write(seq_record, out_handle, fileType)
    qualityCheck(out_file, fileType, in_file, out_file)


def appendNumberToHeaders(in_file):
    #Appends number to each FASTA Header ID
    fileType = fastaOrFastq(in_file)
    with open(in_file, 'r') as in_handle:
        out_file = ''
        for i in in_file.split('.')[0:-1]:
            out_file += i
        out_file += '.unique.' + fileType
        with open(out_file, 'a') as out_handle:
            counter = 0
            for seq_record in SeqIO.parse(in_handle, fileType):
                seq_record.id =  '_' + str(counter)
                SeqIO.write(seq_record, out_handle, fileType)
                counter += 1
    qualityCheck(out_file, fileType, in_file, out_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'appends a string to to all'\
                                     + ' Header IDs in a FASTA/Q file')
    parser.add_argument('input_file',\
                        default = None,\
                        nargs='?',\
                        help = 'FASTA/Q file to modify')
    choices = parser.add_mutually_exclusive_group(required=True)
    choices.add_argument('--string',\
                         default = None,\
                         nargs='?',\
                         help = 'the string to add to each ID the file')
    choices.add_argument('--unique',
                         action='store_true',
                         help='make each ID unique by appending a number')
    parser.add_argument('--version', '-v',\
                        help = 'prints version and exits',\
                        action = 'store_true')
    args = parser.parse_args()
    
    if args.version:
        print(__version__)
        sys.exit(0)
    elif args.input_file == None:
        print(__doc__)
        sys.exit(0)
    elif args.string is not None:
        appendStringToHeaders(args.input_file, args.string)
    else:
        appendNumberToHeaders(args.input_file)

    sys.exit(0)
