#! /usr/bin/env python

from __future__ import print_function

"""Make table summary of Phylosift taxa_summary.txt"""

import argparse
import os
import sys

__author__ = 'Alex Hyer'
__email__ = 'theonehyer@gmail.com'
__license__ = 'GPLv3'
__maintainer__ = 'Alex Hyer'
__status__ = 'Alpha'
__version__ = '0.0.1a1'


def main(args):
    """Run main program

    Args:
        args (NameSpace): argparse NameSpace of variables influencing program
    """

    summary = {}
    args.taxa_summary.readline()
    for line in args.taxa_summary:
        cols = line.strip().split('\t')
        if cols[3] not in summary.keys():
            summary[cols[3]] = {}
        if cols[4] not in summary[cols[3]]:
            summary[cols[3]][cols[4]] = 0
        summary[cols[3]][cols[4]] += float(cols[5])
    for rank in summary.keys():
        rank_total = 0.0
        for key, value in summary[rank].items():
            rank_total += value
        for key, value in summary[rank].items():
            summary[rank][key] = value / rank_total
    for rank in summary.keys():
        for key, value in summary[rank].items():
            args.output.write('{0}\t{1}\t{2}{3}'.format(rank, key,
                                                        str(value),
                                                        os.linesep))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.
                                     RawDescriptionHelpFormatter)
    parser.add_argument('-t', '--taxa_summary',
                        default=sys.stdin,
                        type=argparse.FileType('r'),
                        help='Phylosift taxa_summary.txt file')
    parser.add_argument('-o', '--output',
                        default=sys.stdout,
                        type=argparse.FileType('w'),
                        help='output file for results')
    args = parser.parse_args()

    main(args)

    sys.exit(0)
