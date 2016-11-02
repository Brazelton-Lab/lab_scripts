#! /usr/bin/env python
"""
Convert the results of a database search from one file format to another

Copyright:

    convert_query_format.py Convert query format between databases
    Copyright (C) 2016  William Brazelton

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
from Bio import SearchIO
import textwrap
import re

def io_check(infile, mode='rU'):
    try:
        fh = open(infile)
    except IOError:
        try:
            fh = open(infile, mode)
        except IOError as e:
            print(e)
            sys.exit(1)
    else:
        if mode == 'w':
            print("{} already exists in cwd)".format(os.path.basename(infile)))
            sys.exit(1)
        fh.close()
    return infile

def parse_kwargs(format_arguments):
    r = re.compile("(?P<key>):(?P<value>)")
    keywords = {}
    for argument in format_arguments:
        matched = r.search(argument)
        if not matched:
            print("")
            sys.exit(1)
        kewords[matched.group('key')] = matched.group('value')
    return keywords

def main():
    extensions = {'blast-tab': ['tsv', 'csv', 'blast', 'm8', 'blastm8'],
        'blast-text': ['txt', 'bls', 'blast'], 'blast-xml': ['xml'],
        'blat-psl': ['psl'], 'hmmer3-tab': ['tsv', 'csv'],
        'hmmer3-text': ['txt'], 'hmmer2-text': ['txt'],
        'exonerate-text': ['txt']}
    kwargs = args.keywords
    infile = args.infile
    in_type = args.in_type
    in_ext = infile.split('.')[-1]
    proper_ext = extensions[in_type][0]
    if in_ext not in extensions[in_type]:
        print(textwrap.fill("error: invalid input file extension \"{}\". An "
            "appropriate extension for this input type is {}"
            .format(in_ext, proper_ext), 79))
        sys.exit(1)

    out_type = args.out_type
    if args.output:
        outfile = io_check(args.output, 'w')
    else:
        out_ext = extensions[out_type][0]
        outfile = io_check("{}.{}".format('.'.join(infile.split('.')[:-1]), out_ext), 'w')

    print("output will be in {} and formatted as {}".format(outfile, out_type))
    SearchIO.convert(infile, in_type, outfile, out_type, out_kwargs=kwargs)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="convert the results of a "
                        "database search from one file format to another")
    parser.add_argument('infile', metavar='<input file>',
                        type=io_check,
                        help="input file")
    parser.add_argument('--input-format', metavar='TYPE',
                        dest='in_type',
                        default='blast-text',
                        choices=['blast-text', 'blast-tab', 'blast-xml',
                                 'hmmer3-text', 'hmmer3-tab', 'hmmer2-text',
                                 'exonerate-text', 'blat-psl'],
                        help="input file format. Available options are "
                        "blast-text, blast-tab, blast-xml, blat-psl, "
                        "hmmer3-tab, hmmer3-text, and hmmer2-text [default: "
                        "blast-text]")
    parser.add_argument('--output-format', metavar='TYPE',
                        dest='out_type',
                        default='blast-xml',
                        choices=['blast-tab', 'blast-xml', 'blat-psl',
                                 'hmmer3-tab'],
                        help="output file format. Available options are "
                        "blast-tab, blast-xml, blat-psl, and hmmer3-tab "
                        "[default: blast-xml]")
    parser.add_argument('-k', '--keywords', metavar='KEY',
                        nargs='+',
                        type=parse_kwargs,
                        help="additional format-specific arguments to be "
                        "passed to the converter [ex: -k key:value key:value]")
    parser.add_argument('-o', '--output', metavar='FILE',
                        help="output file name")
    args = parser.parse_args()
    main()
    sys.exit(0)
