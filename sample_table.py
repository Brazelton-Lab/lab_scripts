#! /usr/bin/env python3
"""
Generates new sample tables and appends to a master table based on supplied
inputs.

Copyright:

    sample_tables generates sample tables
    Copyright (C) 2018  William Brazelton

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
from datetime import date
import os
import re
import sys
from time import sleep


def parse_commas(argument):
    """Return list from command line argument with values separated by a comma.
    """
    argument = argument.strip().replace(' ' , '')
    argument = argument.split(',')
    return(argument)


def get_new_id(infile, prefix, sep=','):
    """return the next number in the series from the most recent identifying
    number in the master table."""
    ID = re.compile("{}[0-9]+".format(prefix))

    with open(infile) as in_h:
        header = in_h.readline()
        header = header.split(sep)
        try:
            index = header.index('SampleID')
        except ValueError:
            print("Error: no column named 'SampleID' found in {}"\
                  .format(infile), file=sys.stderr)
            leave = input('Press any key to exit')
            sys.exit(1)

        samples = []
        for line in in_h:
            line = line.strip().split(sep)
            sample = line[index]

            # Verify that the format of the sample ID matches what is expected
            if not ID.match(sample):
                print("Error: sample IDs must consist of one or more "
                      "alphabetical characters followed by one or more "
                      "integers", file=sys.stderr)
                leave = input('Press any key to exit')
                sys.exit(1)

            # Store ID numbers
            try:
                sample_number = int(sample[len(prefix):])
            except ValueError:
                print("Error: sample IDs must consist of one or more "
                      "alphabetical characters followed by one or more "
                      "integers", file=sys.stderr)
                leave = input('Press any key to exit')
                sys.exit(1)
            samples.append(sample_number)

    if len(samples) <= 0:
        new = 1
    else:
        new = sorted(samples)[-1] + 1
    return(new)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-m', '--master',
        default='C:\\Users\\user\Documents\LC2018\LC_2018.sample_table.master.csv',
        help="master sample table [default: C:\\Users\\user\Documents\LC2018\\"
             "LC_2018.sample_table.master.csv]")
    parser.add_argument('-o', '--out_dir',
        default='C:\\Users\\user\Desktop\\',
        help="output directory [default: C:\\Users\\user\Desktop\]")
    parser.add_argument('-i', '--id',
        type=int,
        help="first sample ID number to be assigned to the new table. Ignores "
             "the last sample ID number found in the master table")
    parser.add_argument('-f', '--fields',
        dest='fields',
        default="CruiseID,PrivateID,ContainerType,SampleType,SampleQuantity,\
                 StorageMethod",
        type=parse_commas,
        help="comma-separated list of additional field names to include in "
             "the sample table [default: CruiseID, PrivateID, ContainerType, "
             "SampleType, SampleQuantity, StorageMethod]")
    parser.add_argument('-s', '--sep',
        default=',',
        help="field separator character [default: ,]")
    parser.add_argument('-p', '--pad',
        type=int,
        default=5,
        help="pad the sample ID number with zeros to this width [default: 5]")
    parser.add_argument('-c', '--code',
        type=str,
        default='LC',
        help="project code [default: LC]")
    args = parser.parse_args()

    # Current date in ISO format
    current_date = date.today().isoformat().replace('-', '')

    # Check if master table exists and is not empty
    if not os.path.isfile(args.master):
        with open(args.master, 'w') as master_h:
            master_h.write("SampleID,Owner,Date\n")

    # Obtain input from user
    assignee = input('Enter your name: ')
    # Verify that input does not contain separator character
    if args.sep in assignee:
        print("Error: name '{}' cannot contain the same character '{}' that "
              "is used as the field separator".format(assignee, args.sep))
        leave = input('Press any key to exit')
        sys.exit(1)

    nlabel = input('Enter the desired number of labels: ')
    # Verify that input is of correct type
    try:
        nlabel = int(nlabel)
    except ValueError:
        print("Error: number of labels must be an integer value", \
              file=sys.stderr)
        leave = input('Press any key to exit')
        sys.exit(1)

    # Find the sample ID range
    first_id = get_new_id(args.master, args.code, sep=args.sep)
    last_id = first_id + nlabel
    sample_range = "{}-{}".format(first_id, last_id - 1)

    # Write or append to output files
    output = "{}{}_{}_{}.csv".format(args.out_dir, assignee.replace(" ", "_"), \
             sample_range, current_date)

    print("Outputting table of sample IDs in {}".format(output), \
          file=sys.stdout)
    nfields = len(args.fields)
    ncomma = "," * nfields
    with open(output, 'w') as out_h, open(args.master, 'a') as master_h:
        out_h.write('SampleID,Owner,Date,{}\n'.format(','.join(args.fields)))
        for sid in list(range(first_id, last_id)):
            out_h.write('LC{},{},{}{}\n'.format(str(sid).zfill(args.pad), \
                        assignee, current_date, ncomma))
            master_h.write('LC{},{},{}\n'.format(str(sid).zfill(args.pad), \
                           assignee, current_date))

    print("Identifying number of first sample in series: {}"\
          .format(str(first_id).zfill(args.pad)), file=sys.stdout)

    # Allow time to note ID number before exiting
    sleep(2)
    leave = input('Press any key to exit')

if __name__ == "__main__":
    main()
    sys.exit(0)
