#! /usr/bin/env python
"""
transfers taxonomy from .taxonomy file to each OTU in .shared.transposed file
to make shared.transposed file run transpose-table.py on shared file
"""
from __future__ import print_function
import sys
import os
import argparse
import textwrap

def argument_parser():
    parser = argparse.ArgumentParser(description="add taxonomy information to "
                                     "OTUs in shared file")
    parser.add_argument('taxonomy', metavar='TAX',
                        type=io_check,
                        help="mothur-formatted otu taxonomy file")
    parser.add_argument('shared', metavar='SHARE', 
                        type=io_check,
                        help="transposed shared file")
    return parser

def io_check(infile, mode='rU'):
    try:
        fh = open(infile, mode)
        fh.close()
    except IOError as e:
        print(e)
        sys.exit(1)
    return infile

args = argument_parser().parse_args()
outfile = io_check("{}.taxonomy".format(args.shared), 'w')

D = {}
with open(args.taxonomy) as taxfile:
    header = taxfile.readline()
    num_cols = len(header.split('\t'))
    if num_cols != int(3):
        print(textwrap.fill("error: {}: incorrect file format. The OTU taxonomy file should have three columns. Actual number of columns: {}".format(os.path.basename(args.taxonomy), str(num_cols)), 79))
        sys.exit(1)
    if not header.lower().startswith("otu"):
        taxfile.seek(0)
    for line in taxfile:
        cols = line.strip().split('\t')
        otu = cols[0]
        tax = cols[-1].replace(';', '\t')
        D[otu] = tax

with open(outfile, 'w') as out:
    with open(args.shared) as sharedfile:
        for line in sharedfile:
            if "label" in line:
                out.write(line)
            elif "Group\t" in line:
                num_col = len(line.strip().split())
                out.write(line)
            elif "numOtus\t" in line: 
                out.write(line)
            else:
                line = line.strip().split()
                if len(line) != num_col:
                    continue
                shared_otu = line[0]
                try: 
                    tax = D[shared_otu]
                except KeyError as e:
                    print("error: can't find {} in {}".format(e, os.path.basename(args.shared)))
                    sys.exit(1)
                out.write("{}\t{}\n".format('\t'.join(line), tax))
