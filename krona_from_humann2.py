#!/usr/bin/env python
"""
Rearranges HUMAnN2 abundance output into a TSV for use with ktImportText

Copyright:

    krona_from_humann2  Rearranges HUMAnN2 abundance output into a TSV for use with ktImportText

    Copyright (C) 2016  William Brazelton <comma-separated list of authors>

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
