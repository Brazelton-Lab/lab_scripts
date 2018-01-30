#! /usr/bin/env python
"""
Create a multiple sample table of Pathway Abundances using individual sample 
pathabundance files from HUMAnN2

Copyright:

    table_from_humann2  Create a multiple sample table of Pathway Abundances using individual sample path abundance files from HUMAnN2

    Copyright (C) 2016  William Brazelton

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.﻿
"""

from __future__ import print_function
import argparse
import sys
import os

def print_path(pathway, database, size):
    abundances = database[pathway]
    diff = size - len(abundances)
    if diff > 0:
        abundances.extend(list('0') * diff)
    print("{}\t{}".format(pathway, '\t'.join(abundances)))

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('infiles', metavar='pathabundance',
        nargs='+',
        help="space-separated list of HUMAnN2 pathabundance files")

    parser.add_argument('-l', '--labels', metavar='"label, label, ..."',
        type=str,
        help="sample names to use in the report, given in a comma-separated "
            "list [default: use file names]. The order given must match the "
            "order in which the files are provided. If spaces are used "
            "in-between labels, the full argument must be wrapped in quotes.")

    args = parser.parse_args()

    labels = args.labels.split(',')
    labels = [i.lstrip() for i in labels]

    valid_files = []
    # check that the can be opened
    for index, infile in enumerate(args.infiles):
        infile = args.infiles[index]
        try:
            fh = open(infile).close()
        except IOError:
            basename = os.path.basename(infile)
            print("Unable to open {} ... skipping".format(basename), file=sys.stderr)
            if labels:
                labels = labels[:index] + labels[index + 1:]
        else:
            valid_files.append(infile)

    path_db = {}

    for index, infile in enumerate(valid_files):
        with open(infile) as in_h:
            for line in in_h:
                if line.startswith('#'):
                    continue
                try:
                    pathway, abundance = line.strip().split('\t')
                except ValueError:
                    basename = os.path.basename(infile)
                    print("file '{}' does not have a recognizable format".format(basename), file=sys.stderr)
                    sys.exit(1)
                if pathway in path_db:
                    diff = index - len(path_db[pathway])
                    if diff > 0:
                        path_db[pathway].extend(list('0') * diff)
                    path_db[pathway].append(abundance)
                else:
                    path_db[pathway] = list('0') * index + [abundance]

    num_valid = len(valid_files)

    if labels:
        header = "# Pathway\t{}".format('\t'.join(labels))
    else:
        header = "# Pathway\t{}".format('\t'.join(valid_files))
    print(header)

    if 'UNMAPPED' in path_db:
        print_path('UNMAPPED', path_db, num_valid)

    if 'UNINTEGRATED' in path_db:
        print_path('UNINTEGRATED', path_db, num_valid)

    for pathway in sorted(path_db):
        if pathway == 'UNMAPPED' or pathway == 'UNINTEGRATED':
            continue
        print_path(pathway, path_db, num_valid)

if __name__ == "__main__":
    main()
    sys.exit(0)
