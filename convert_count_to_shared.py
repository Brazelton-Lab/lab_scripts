#!/usr/bin/env python

from __future__ import print_function

"""Converts a MOTHUR count file to a MOTHUR shared file

Usage:
    convert_count_to_shared [count_file]

Synopsis:
    The generation of OTUs are not always desired but many functions performed
    by MOTHUR require the OTU generation step. For users wishing to take
    advantage of these functions without OTU generation,
    convert_count_to_shared creates a '.shared' file from a '.count_table'
    file where each Representative Sequence is labeled as a unique OTU. This
    effectively bypasses OTU generation.

Required Arguments:
    count_file:     MOTHUR generated count table file to convert to a shared
                    file

Output:
    A MOTHUR formatted '[count_file].shared' file.
"""

__version__ = '1.0.0'
__author__ = 'Alex Hyer, William Brazelton'

import argparse
from bio_utils.file_tools.file_check import FileChecker
import sys


def read_count_file(count_file):
    with open(count_file, 'rU') as count_handle:
        first_line = count_handle.readline()
        first_line = first_line.strip().split('\t')
        count_table = {header: [] for header in first_line}
        for line in count_handle:
            split_line = line.strip().split('\t')
            for column in enumerate(split_line):
                header = first_line[column[0]]
                count_table[header].append(column[1])
    return count_table


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.
                                     RawDescriptionHelpFormatter)
    parser.add_argument('count_file', metavar='Count File',
                        help='MOTHUR count file to convert to shared file')
    args = parser.parse_args()

    count_file = open(args.count_file)
    out_name = count_file.name + '.shared'

    count_table = read_count_file(count_file.name)
    with open(out_name.name, 'w') as out_handle:
        new_first_line = 'label\tGroup\tnumOtus\t' + '\t'.join(
            [sequence for sequence in count_table['Representative_Sequence']])
        out_handle.write(new_first_line + '\n')
        for key in count_table.keys():
            if key == 'Representative_Sequence' or key == 'total':
                continue
            column_number = len(count_table[key])
            new_line = 'unique\t{0}\t{1}\t'.format(key, column_number) \
                       + '\t'.join(abundance for abundance in count_table[key])
            out_handle.write(new_line + '\n')

if __name__ == '__main__':
    main()
sys.exit(0)
