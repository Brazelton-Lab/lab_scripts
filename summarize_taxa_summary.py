#! /usr/bin/env python

"""
Make table summary of Phylosift taxa_summary.txt

Copyright:

    summarize_taxa_summary  Make table summary of Phylosift taxa_summary.txt

    Copyright (C) 2016  William Brazelton, Alex Hyer

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from __future__ import print_function
import argparse
import os
import sys

__author__ = 'Alex Hyer'
__email__ = 'theonehyer@gmail.com'
__license__ = 'GPLv3'
__maintainer__ = 'Alex Hyer'
__status__ = 'Production'
__version__ = '1.0.0'


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
            summary[rank][key] = value / rank_total * 100.0
    for rank in summary.keys():
        if args.greatest is False:
            for key, value in summary[rank].items():
                args.output.write('{0}\t{1}\t{2}{3}'.format(rank, key,
                                                            str(value),
                                                            os.linesep))
        else:
            greatest = max(summary[rank].items(), key=lambda x: x[1])
            args.output.write('{0}\t{1}\t{2}{3}'.format(rank, greatest[0],
                                                        str(greatest[1]),
                                                        os.linesep))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.
                                     RawDescriptionHelpFormatter)
    parser.add_argument('-g', '--greatest',
                        action='store_true',
                        help='only report greatest value')
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
