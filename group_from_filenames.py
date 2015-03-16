#! /usr/bin/env python

from __future__ import print_function

import sys
import argparse

def argument_parser():
    parser = argparse.ArgumentParser(description="create a group file from a "
                                     "text file containing a newline-separated "
                                     "list of fasta file names. The fasta "
                                     "files should be located in the same "
                                     "directory as the text file and the titles "
                                     "should contain the name of the sample"
                                     "")
    parser.add_argument('infile', metavar='FILE',
                        help="file containing a newline-separated list of "
                        "fasta files with sample names in the title")
    parser.add_argument('-s', '--separator', metavar='CHAR',
                        type=str,
                        default='.',
                        help="specifiy the character that spearates the sample "
                        "name from the rest of the file name [default .]")
    parser.add_argument('-p', '--position', metavar='POS',
                        type=int,
                        default=1,
                        help="specify the location of the sample name in the "
                        "file names [default: 1]")
    return parser

def file_check(infile, mode):
    try:
        fh = open(infile, mode)
        fh.close()
    except IOError as e:
        print(e)
        sys.exit(1)

def parse_file(infile):
    files = []
    file_check(infile, 'rU')
    with open(infile, 'rU') as in_h:
        for line in in_h:
            sample_file = line.strip()
            file_check(sample_file, 'rU')
            files.append(sample_file)
    return files

def main():
    args = argument_parser().parse_args()
    sample_files = parse_file(args.infile)
    outfile = 'all.groups'
    file_check(outfile, 'w')
    with open(outfile, 'w') as out:
        for sample_file in sorted(sample_files):
            sample = sample_file.split(args.separator)[args.position - 1]
            with open(sample_file, 'rU') as s_h:
                for line in s_h:
                    if not line.startswith('>'):
                        continue
                    ident = line.strip('>\n')
                    if ':' in ident:
                        ident = ident.replace(':', '_')
                    output = "{}\t{}\n".format(ident, sample)
                    out.write(output)

if __name__ == "__main__":
    main()
