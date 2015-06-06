#! /usr/bin/env python
"""
Convert the results of a database search from one file format to another
"""

from __future__ import print_function
import sys
import os
import argparse
from Bio import SearchIO
import textwrap

def io_check(infile, mode='rU'):
    try:
        fh = open(infile)
    except IOError:
        try:
            fh = open(infile, mode)
        except IOErroras e:
            print(e)
            sys.exit(1)
    else:
        if mode == 'w':
            print("{} already exists in cwd)".format(os.path.basename(infile)))
            sys.exit(1)
        fh.close()
    return infile

def parse_kwargs(format_arguments):
    keywords = {}
    
    return keywords

def main():
    exts = {'blast-tab': 'csv', 'blast-text': 'txt', 'blast-xml': 'xml',
            'blat-psl': 'psl', 'hmmer3-tab': 'csv', 'hmmer3-text': 'txt', 
            'hmmer2': 'txt', 'exonerate-text': 'txt'}
    kwargs = args.options
    infile = args.infile
    in_type = args.in_type
    in_ext = infile.split('.')[-1]
    proper_ext = exts[in_type]
    if in_ext != exts[in_type]:
        print(textwrap.fill("error: invalid input file extension \"{}\". The "
            "appropriate extension for this input type is {}"
            .format(in_ext, proper_ext), 79))
        sys.exit(1)
    out_type = args.out_type
    out_ext = ext[out_type]
    outfile = io_check("{}.{}".format('.'.join(infile.split('.')[:-1]), out_ext))
    SearchIO.convert(infile, in_type, outfile, out_type, out_kwargs=kwargs)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="convert the results of a "
                        "database search from one file format to another")
    parser.add_argument('infile', metavar='<input blast>',
                        type=io_check,
                        help="")
    parser.add_argument('--in-type', metavar='TYPE',
                        dest='in_type',
                        default='blast-text'
                        choices=['blast-text','blast-tab','blast-xml','hmmer3-text','hmmer3-tab','hmmer2-text','exonerate-text','blat-psl'],
                        help="input file format [default: blast-text]")
    parser.add_argument('--out-type', metavar='TYPE',
                        dest='out_type',
                        default='blast-xml',
                        choices=['blast-tab','blast-xml','blat-psl','hmmer3-tab'],
                        help="output file format [default: blast-xml]")
    parser.add_argument('-k', '--keywords',
                        nargs='+',
                        type=parse_kwargs,
                        help="[ex: --op")
    args = parser.parse_args()
    main()
    sys.exit(0)
