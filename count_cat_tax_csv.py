#! /usr/bin/env python

"""Adds taxonomy information from CSV files produced by phyloseq package (for example) to otu tables

Copyright:

    count_cat_tax_csv.py Combine taxonomy file and count table
    Copyright (C) 2017  William Brazelton, Alex Hyer

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

import argparse
import sys


def main():
    otus = {}
    header = [',Domain,Phylum,Class,Order,Family,Genus']
    with open(args.count_table, 'rU') as count_handle:
        header = header + (count_handle.readline().strip().lstrip(',').split(','))
        for line in count_handle:
            columns = line.strip().split(',')
            otus[columns[0]] = [','.join(columns[1:])]
    with open(args.tax_table, 'rU') as tax_handle:
        tax_handle.readline().strip().split(',')
	for line in tax_handle:
            columns = line.strip().split(',')
            otus[columns[0]] = columns[1:] + otus[columns[0]]
    with open(args.output, 'w') as out_handle:
	for i in header:        
		out_handle.write(str(i) + ',')
	out_handle.write('\n')
        for key in otus.keys():
            if len(otus[key]) == 1:
                print('WARNING: {0} may not be in {1}'.format(
                      key, args.count_table))
                otus[key] = ['NA' for i in range(len(header[1:]))] + otus[key]
            output = key + ',' + ','.join(otus[key]) + '\n'
            out_handle.write(output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                    formatter_class=argparse.
                                    RawDescriptionHelpFormatter)
    parser.add_argument('-t', '--tax_table',
                        required=True,
                        help='taxonomy CSV file from e.g. phyloseq')
    parser.add_argument('-c', '--count_table',
                        required=True,
                        help='OTU / count table file from e.g. phyloseq')
    parser.add_argument('-o', '--output',
                        required=True,
                        help='File to write new table to')
    args = parser.parse_args()
    main()
    sys.exit(0)
