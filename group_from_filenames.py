#! /usr/bin/env python

from __future__ import print_function

import sys
import os
import argparse
import textwrap

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
            for line in in_h:
                if not line.startswith('>'):
                    continue
                identifier = line.strip('>\n').split()[0]
                if ':' in identifier:
                    identifier = identifier.replace(':', '_')
                print("{}\t{}".format(identifier, sample))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create a group file for a "
        "project by providing one or more fasta files. The file names should "
        "contain the name of the sample.")
    parser.add_argument('infiles', metavar='fasta_file',
        nargs='+',
        type=io_check,
        help="newline-separated list of fasta files")
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
