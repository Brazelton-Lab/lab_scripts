#! /usr/bin/env python
# sort sequences in standard.fasta according to bar code
# files saved with filenames according to barcodes_distribution.txt

from __future__ import print_function
import sys
import os
import re
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="creates a mothur-formatted "
        "group file from a fasta file where the sample names are in the "
        "sequence ids")
    parser.add_argument('fasta',
                        type=io_check,
                        help="input fasta file")
    args = parser.parse_args()
    infile = args.fasta
    outfile = os.path.basename(infile) + '.group'

    r = re.compile("(?<=[Ss][Ee][Rr][Pp]_)(?P<name>SE_\w+)(?=\|)")
    with open(infile, 'rU') as in_h:
        with open(outfile, 'w') as out_h:
            for line in in_h:
                if line.startswith('>'):
                    header = line.strip('>\n')
                    header = header.split()[0]
                    matched = r.search(header)
                    if not matched:
                        print("error: sequence ids are not formatted in such a way that this program can parse")
                        sys.exit(1)
                    sample = matched.group('name')
                    out_h.write("{}\t{}\n".format(header, sample))
                else:
                    continue 
    sys.exit(0)
