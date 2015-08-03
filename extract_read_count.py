#! /usr/bin/env python
"""
For extracting read count information from idba_ud assemblies
"""

from __future__ import print_function
import sys
import argparse

def main():
    infile = args.infile
    try:
        fh = open(infile, 'rU')
    except IOError as e:
        print(e)
        sys.exit(1)
    else:
        fh.close()

    with open(infile) as in_h:
        for line in in_h:
            if line.startswith('>'):
                header = line.strip('>\n').split(' ')
                contig_id = header[0]
                read_count = header[-1].split('_')[-1]
                print("{}\t{}".format(contig_id, read_count))
    sys.exit(0)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="pull out read count information from an idba_ud assembly "
        "into a tab-delimited file. Output is to STDOUT. To put into a file, "
        "use redirection (ex: extract_read_count.py <contigs.fasta> > "
        "<outfile>")
    parser.add_argument('infile', metavar="contigs", 
        help="contig file from idba_ud")
    args = parser.parse_args()
    main()
