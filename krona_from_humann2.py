#!/usr/bin/env python
"""Rearranges HUMAnN2 abundance output into a TSV for use with ktImportText"""
from __future__ import print_function

import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('abund_file', metavar='abundance_file',
        help='HUMAnN2 pathways abundance file in tsv format')

    args = parser.parse_args()

    with open(args.abund_file) as in_h:
        for line in in_h:
            # Skip header
            if line.startswith('#'):
                continue

            try:
                pathways, abundance = line.strip().split('\t')
            except ValueError:
                print("Unknown input file format", file=sys.stderr)
                sys.exit(1)

            pathways = pathways.split(';')
            if pathways != ['UNMAPPED'] and pathways != ['UNINTEGRATED']:
                pathways.insert(0, 'INTEGRATED')

            print('{!s}\t{}\n'.format(abundance, '\t'.join(pathways)))

if __name__ == '__main__':
    main()
    sys.exit(0)
