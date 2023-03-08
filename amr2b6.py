#! /usr/bin/python3.6
"""
Convert the output of AMRFinder to B6 format.

Required input is the tab-separated best-hit results of AMRFinder searches to 
NCBIfam-AMR databases. Optional input is a relational database in JSON format,  
which can be used to replace the subject identifier from the results to a new 
identifier for the purpose of annotation.

The compression algorithm is automatically detected for input files based on
the file extension. To compress output, add the appropriate file extension to
the output file name (e.g. .gz, .bz2). Leave off '--out' to direct output to
standard output (stdout). Standard input (stdin) can be redirected to the 
input file if supplied with '-'.

Copyright:

    amr2b6 Convert AMRFInder output to B6 format
    Copyright (C) 2022  William Brazelton

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
from arandomness.argparse import Open, ParseSeparator
from seq_annot import seqio, db
import sys
import textwrap
from time import time

__author__ = "Christopher Thornton"
__license__ = 'GPLv3'
__maintainer__ = 'Christopher Thornton'
__status__ = "Alpha"
__version__ = "0.0.2"


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('csv',
        metavar='in.csv',
        help="input tabular results of database searches from AMRFinder")
    parser.add_argument('-m', '--mapping',
        metavar='in.json',
        dest='map_file',
        type=str,
        help="input relational database in JSON format")
    parser.add_argument('-f', '--field',
        type=str,
        help="field in the relational database corresponding to the old ID")
    parser.add_argument('-o', '--out',
        metavar='out.b6',
        action=Open,
        mode='wb',
        default=sys.stdout,
        help="output results in B6 format")
    parser.add_argument('--split',
        action='store_true',
        help="split field value from relational database entry on first "
            "encounter of ':' for compatibility with certain feature "
            "attributes such as Dbxref. The second item of the split value "
            "will be stored as the new ID")
    args = parser.parse_args()

    # Input sanity check
    if (args.map_file and not args.field) or (args.field and not args.map_file):
        parser.error("arguments -m/--mapping and -f/--field must be supplied "
            "together")

    # Output run information
    all_args = sys.argv[1:]
    print("{} {!s}".format('amr2b6', __version__), file=sys.stderr)
    print(textwrap.fill("Command line parameters: {}"\
          .format(' '.join(all_args)), 79), file=sys.stderr)

    # Track program run-time
    start_time = time()

    # Assign variables based on user inputs
    if args.csv == '-':
        in_h = seqio.open_io(sys.stdin, mode='rb')
    else:
        in_h = seqio.open_io(args.csv, mode='rb')

    out_h = args.out

    field = args.field

    # Load mappings, if provided
    mapping = db.load_dbs([args.map_file]) if args.map_file else {}

    # Invert Keys and Values
    map_invert = {}
    for entry_id in mapping:
        try:
            new_id = mapping[entry_id][field]
        except KeyError:
            print("warning: field {} not found for entry {} "\
                .format(field, entry_id), file=sys.stderr)
            continue

        if args.split:
            try:
                final_id = new_id.split(':', 1)[1]
            except ValueError:
                final_id = new_id
        else:
            final_id = new_id

        if not final_id in map_invert:
            map_invert[final_id] = entry_id
        else:
            print("warning: relationship between field value {} and entry ID "
                "not one-to-one for entry {}. Entry will be skipped"\
                .format(new_id, entry_id), file=sys.stderr)
    
    # Convert output to B6
    header = in_h.readline()
    header = header.decode('utf-8')
    header = header.rstrip().split('\t')
    header_map = {j: i for i, j in enumerate(header)}

    hit_totals = 0
    for nline, line in enumerate(in_h):
        line = line.decode('utf-8')

        if line.startswith('#'):  #skip comments
            continue
        else:
            hit_totals += 1

        line = line.rstrip().split('\t')

        target = line[header_map['Protein identifier']]
        subject = line[header_map['HMM id']]
        pident = line[header_map['% Identity to reference sequence']] if \
            line[header_map['% Identity to reference sequence']] != "NA" else "-"
        alength = line[header_map['Alignment length']] if \
            line[header_map['Alignment length']] != "NA" else "-"

        # Convert subject ID if mapping file provided
        if mapping:
            try:
                subject_id = map_invert[subject]
            except KeyError:
                print("error: {}, line {}: subject {} does not have a "
                    "corresponding entry for field {} in database {}"\
                    .format(args.csv, nline, subject, field, args.map_file), \
                    file=sys.stderr)
                print(map_invert, file=sys.stderr)
                sys.exit(1)
        else:
            subject_id = subject

        seqio.write_io(out_h, "{}\t{}\t{}\t{}\t-\t-\t-\t-\t-\t-\t-\t-\n"\
            .format(target, subject_id, pident, alength))

    # Calculate and print program run-time
    end_time = time()
    total_time = (end_time - start_time) / 60.0
    print("", file=sys.stderr)
    print("It took {:.2e} minutes to convert {!s} hits".format(total_time, \
        hit_totals), file=sys.stderr)
    print("", file=sys.stderr)


if __name__ == "__main__":
    main()
    sys.exit(0)
