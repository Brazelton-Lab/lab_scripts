#! /usr/bin/env python

"""

Copyright:

    group_from_filenames  Creates a group file for a project by providing one or more fasta or fastq formatted files.

    Copyright (C) 2016  William Brazelton <comma-separated list of authors>

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

from __future__ import print_function

import sys
import os
import argparse
import textwrap
import screed

def io_check(infile):
    try:
        fh = open(infile, 'rU')
    except IOError as e:
        print(textwrap.fill("Error: unable to read {}".format(
            os.path.basename(infile)), 79))
        sys.exit(1)
    else:
        fh.close()
    return infile

def print_out(header, sample):
    print("{}\t{}".format(header.replace(':', '_').split(' ')[0], sample))

def main():
    for infile in sorted(args.infiles):
        file_name = os.path.basename(infile).split(args.separator)
        parts = len(file_name)
        if args.position > parts + 1:
            print(textwrap.fill("Error: invalid position \"{}\". There are "
                "only {} groups created when splitting {} by {}".format(
                str(args.position), str(parts), os.path.basename(infile), 
                args.separator), 79))
            sys.exit(1)
        sample = file_name[args.position - 1]
        with open(infile, 'rU') as in_h:
            line = in_h.readline()
            in_h.seek(0)
            if line.startswith('>'):
                for seq in screed.fasta.fasta_iter(in_h):
                    print_out(seq.name, sample)
            elif line.startswith('@'):
                for seq in screed.fastq.fastq_iter(in_h):
                    print_out(seq.name, sample)
            else:
                print("Error: unrecognized file format")
                sys.exit(1)
    sys.exit(0)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create a group file for a "
        "project by providing one or more fasta or fastq formatted files. The "
        "file names should contain the name of the sample.")
    parser.add_argument('infiles', metavar='seq_file',
        nargs='+',
        type=io_check,
        help="space-separated list of fasta or fastq files")
    parser.add_argument('-s', '--separator', metavar='CHAR',
        default='.',
        help="character that separates the sample name from the rest of the "
        "file name [default .]")
    parser.add_argument('-p', '--position', metavar='INT',
        type=int,
        default=1,
        help="location of the sample name in the file name [default: 1]")
    args = parser.parse_args()
    main()
