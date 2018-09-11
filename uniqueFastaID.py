#!/usr/bin/env python

"""uniqueFastaID alters the ID of each FASTA entry in a file

Some programs such as IDBA-UD may create multiple fasta/q files where
the ID of the entries in each file are identical. uniqueFastaId simply takes a
fasta/q file and appends a user defined string to each Fasta Header ID to
ensure that each file will not have headers that clash if the files need to be
concatenated for use by various programs. Unqique adds a number to each ID to
avoid repeats.

Copyright:

    uniqueFastaId  Makes the id of each FASTA entry in a file unique

    Copyright (C) 2016  William Brazelton, Alex Hyer

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import argparse
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from bioinformatic_tools import fastaOrFastq
from bioinformatic_tools import qualityCheck

__version__ = '0.12'


def appendStringToHeaders(in_file, string_to_append):
    # Appends string to each FASTA Header ID
    fileType = fastaOrFastq(in_file)
    with open(in_file, 'r') as in_handle:
        out_file = ''
        for i in in_file.split('.')[0:-1]:
            out_file += '{0}.'.format(i)
        out_file += 'string.' + in_file.split('.')[-1]
        with open(out_file, 'a') as out_handle:
            for seq_record in SeqIO.parse(in_handle, fileType):
                if args.replace is True:
                    seq_record.id = ''
                if args.strip_description is True:
                    seq_record.description = ''
                seq_record.id += string_to_append
                SeqIO.write(seq_record, out_handle, fileType)
    qualityCheck(out_file, fileType, in_file, out_file)


def appendNumberToHeaders(in_file):
    # Appends number to each FASTA Header ID
    fileType = fastaOrFastq(in_file)
    with open(in_file, 'r') as in_handle:
        out_file = ''
        for i in in_file.split('.')[0:-1]:
            out_file += '{0}.'.format(i)
        out_file += 'unique.' + in_file.split('.')[-1]
        with open(out_file, 'a') as out_handle:
            counter = 0
            for seq_record in SeqIO.parse(in_handle, fileType):
                if args.replace is True:
                    seq_record.id = ''
                seq_record.id += str(counter)
                if args.strip_description is True:
                    seq_record.description = ''
                SeqIO.write(seq_record, out_handle, fileType)
                counter += 1
    qualityCheck(out_file, fileType, in_file, out_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='appends a string to to all '
                                     'Header IDs in a FASTA/Q file')
    parser.add_argument('input_file',
                        default=None,
                        nargs='*',
                        help='FASTA/Q file to modify')
    parser.add_argument('--replace',
                        action='store_true',
                        help='replace header ID with string or number')
    parser.add_argument('--strip_description',
                        action='store_true',
                        help='remove FASTA header description')
    choices = parser.add_mutually_exclusive_group(required=True)
    choices.add_argument('--string',
                         default=None,
                         nargs='?',
                         help='the string to add to each ID the file')
    choices.add_argument('--unique',
                         action='store_true',
                         help='make each ID unique by appending a number')
    parser.add_argument('--version', '-v',
                        help='prints version and exits',
                        action='store_true')
    args = parser.parse_args()
    
    if args.version:
        print(__version__)
        sys.exit(0)
    elif args.input_file is None:
        print(__doc__)
        sys.exit(0)
    elif args.string is not None:
        for f in args.input_file:
            appendStringToHeaders(f, args.string)
    else:
        for f in args.input_file:
            appendNumberToHeaders(f)

    sys.exit(0)
