#!/usr/bin/env python

from __future__ import print_function

"""Extend elements of a table w/ s given string to match the
longest line of the table"""

import argparse
import os
from StringIO import StringIO
import sys


def max_line(file_handle, delimiter):
    longest_line = 0
    for line in file_handle:
        line_length = len(line.split(delimiter))
        if line_length > longest_line:
            longest_line = line_length
    file_handle.seek(0)
    return longest_line


def smooth_table(in_handle, out_handle, modifier, delimiter):
    longest_line = max_line(in_handle, delimiter)
    for line in in_handle:
        elements = line.strip().split(delimiter)
        while len(elements) < longest_line:
            elements.append(modifier)
        elements[-1] = elements[-1] + os.linesep
        out_handle.write(delimiter.join(elements))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.
                                     RawDescriptionHelpFormatter)
    parser.add_argument('--in_table', '-i',
                        default=sys.stdin,
                        help='Table to modify [Default: STDIN]')
    parser.add_argument('--out_table', '-o',
                        default=sys.stdout,
                        help='File to write modified table to '
                             '[Default: STDOUT]')
    parser.add_argument('--modifier', '-m',
                        default='UNDEFINED',
                        help='String to place in added table elements '
                             '[Default: UNDEFINED]')
    parser.add_argument('--delimiter', '-d',
                        default='\t',
                        help='Character delimiting table [Default: \\t]')
    args = parser.parse_args()

    try:
        with open(args.in_table, 'rU') as in_handle, \
             open(args.out_table, 'w') as out_handle:
            smooth_table(in_handle, out_handle, args.modifier, args.delimiter)
    except TypeError:
        std_input = ''
        for line in sys.stdin:
            std_input += line
        args.in_table = StringIO(std_input)
        smooth_table(args.in_table, args.out_table, args.modifier,
                     args.delimiter)         


if __name__ == '__main__':
    main()
    sys.exit(0)
