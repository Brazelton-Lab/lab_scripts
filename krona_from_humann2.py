#!/usr/local/env python

"""Rearranges HumanN2 abundance output into a TSV for use with ktImportText"""

import argparse
import os
import sys

def main():
    outFile = os.getcwd() + os.sep + args.abundance_file + '.krona'
    with open(args.abundance_file, 'rU') as in_handle:
        with open(outFile, 'w') as out_handle:
            # Skip header
            in_handle.readline()
            for line in in_handle:
                columns = line.strip().split('\t')
                pathways = columns[0].split(';')
                out_handle.write('{0}\t{1}\n'.format(columns[1],
                                 '\t'.join(pathways)))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.
                                     RawDescriptionHelpFormatter)
    parser.add_argument('abundance_file', metavar='HumanN2 Abundance File',
                        help='HumanN2 Abundance TSV File')
    args = parser.parse_args()

    main()
    sys.exit(0)
