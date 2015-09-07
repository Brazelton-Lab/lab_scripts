#! /usr/bin/env python
# sort sequences in standard.fasta according to bar code
# files saved with filenames according to barcodes_distribution.txt

from __future__ import print_function
import sys
import os
import argparse

def io_check(infile, mode='rU'):
    try:
        fh = open(infile, mode)
    except IOError as e:
        print(e)
        sys.exit(1)
    else:
        fh.close()
    return infile

def main():
    parser = argparse.ArgumentParser(description="creates a mothur-formatted "
        "group file from a fasta file where the sample names are in the "
        "sequence ids")
    parser.add_argument('fasta',
                        type=io_check,
                        help="input fasta file")
    parser.add_argument('-s', '--separator', metavar='CHAR',
        default='|',
        help="character that separates the sample name from the rest of the "
        "header [default: | ]")
    parser.add_argument('-p', '--position', metavar='INT',
        type=int,
        default=2,
        help="position of the sample name in the header when separated by "
        "the character separator (starting from 0) [default: 2]")

    args = parser.parse_args()
    infile = args.fasta
    outfile = os.path.basename(infile) + '.group'
    position = args.position
    separator = args.separator

    with open(infile, 'rU') as in_h:
        with open(outfile, 'w') as out_h:
            for line in in_h:
                if line.startswith('>'):
                    header = line.strip('>\n')
                    header = header.split()[0]
                    group = header.split(separator)
                    size = len(group)
                    if size <= position:
                        print("error: the position chosen is outside of range")
                        sys.exit(1)
                    group = group[position]
                    out_h.write("{}\t{}\n".format(header, group))

if __name__ == "__main__":
    main()
    sys.exit(0)
