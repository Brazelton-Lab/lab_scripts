#! /usr/bin/env python

"""Adds taxonomy information from TSV files to count tables
"""

import argparse
import sys


def main():
    otus = {}
    header = ['taxonomy']
    with open(args.tsv, 'rU') as tsv_handle:
        for line in tsv_handle:
            columns = line.strip().split('\t')         
            otus[columns[0]] = [';'.join(columns[1:])]
    with open(args.count_table, 'rU') as count_handle:
        header = count_handle.readline().strip().split('\t') + header
        for line in count_handle:
            columns = line.strip().split('\t')         
            otus[columns[0]] = columns[1:] + otus[columns[0]]
    with open(args.output, 'w') as out_handle:
        out_handle.write('\t'.join(header) + '\n')
        for key in otus.keys():
            output = key + '\t' + '\t'.join(otus[key]) + '\n'
            out_handle.write(output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                    formatter_class=argparse.
                                    RawDescriptionHelpFormatter)
    parser.add_argument('-t', '--tsv',
                        required=True,
                        help='taxonomy TSV file from taxonomy2tsv.py')
    parser.add_argument('-c', '--count_table',
                        required=True,
                        help='Count table file from Mothur')
    parser.add_argument('-o', '--output',
                        required=True,
                        help='File to write new table to')
    args = parser.parse_args()
    main()
    sys.exit(0)
