#! /usr/bin/env python3
"""
Generates sample table and new labels based on the number of desired labels
and the last sample identifier used.

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
import subprocess
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
        default='C:\\Users\\user\\Documents\\LC2018\\LC_2018.sample_table.master.csv',
        help="master sample table [default: C:\\Users\\user\Documents\LC2018\\"
             "LC_2018.sample_table.master.csv]")
    parser.add_argument('-o', '--out-dir',
        metavar='DIR',
        dest='out_dir',
        default='C:\\Users\\user\Desktop\\',
        help="output directory [default: C:\\Users\\user\\Desktop\]")
    parser.add_argument('-i', '--id',
        type=int,
        help="first sample ID number to be assigned to the new table. Ignores "
             "the last sample ID number found in the master table")
    parser.add_argument('-f', '--fields',
        dest='fields',
        default="Date,CruiseID,PrivateID,ContainerType,SampleType,SampleQuantity,\
                 StorageMethod",
        type=parse_commas,
        help="comma-separated list of additional field names to include in "
             "the sample table [default: CruiseID, PrivateID, ContainerType, "
             "SampleType, SampleQuantity, StorageMethod]")
    parser.add_argument('-s', '--sep',
        metavar='CHAR',
        default=',',
        help="field separator character [default: ,]")
    parser.add_argument('-p', '--pad',
        metavar='INT',
        type=int,
        default=5,
        help="pad the sample ID number with zeros to this width [default: 5]")
    parser.add_argument('-c', '--code',
        type=str,
        default='LC',
        help="project code [default: LC]")
    parser.add_argument('-z', '--zpl-dir',
        metavar='DIR',
        dest='zpl_dir',
        default='C:\\Users\\user\\Documents\\ZPL\\',
        help="directory to store the zpl files [default: C:\\Users\\user\\Documents\\ZPL\\]")
    parser.add_argument('-n', '--net',
        metavar='BAT',
        default='C:\\Users\\user\\Documents\\netprint.bat',
        help="location of network printer batch file [default: C:\\Users\\user\\Documents\\netprint.bat]")
    args = parser.parse_args()

    sep = args.sep

    zpl_template = ("^XA\n^FO30,40^ATN^FDLost City 2018^FS\n^FO280,35^FDMM,A{0}{1}"
                   "^BQN,2,6,H^FS\n^FO278,190^ATN^FD{0}{1}^FS\n^FO30,120^ARN^FD{2}"
                   "^FS\n^FO30,172^ARN2^FB250,2,,L^FD{3}^FS\n^XZ\n")

    # Current date in ISO format
    current_date = date.today().isoformat().replace('-', '')

    # verify that network printer batch file exists
    if not os.path.exists(args.net):
        print("Error: cannot find location of the network printer batch file. "
              "Please contact an administrator for assistance")
        leave = input('Press any key to exit')
        sys.exit(1)

    # make output directories if they don't already exist
    for directory in [args.out_dir, args.zpl_dir]:
        if not os.path.exists(directory):
            os.mkdirs(directory)

    # Check if master table exists and is not empty
    if not os.path.isfile(args.master):
        with open(args.master, 'w') as master_h:
            master_h.write("SampleID{0}Owner{0}{1}\n".format(sep, \
                           sep.join(args.fields)))

    # Obtain input from user
    assignee = input('Enter your name, given followed by family (e.g. Jane Doe): ')
    # Verify that input does not contain separator character
    if args.sep in assignee:
        print("Error: name '{}' cannot contain the same character '{}' that "
              "is used as the field separator".format(assignee, sep))
        leave = input('Press any key to exit')
        sys.exit(1)

    owner = assignee.split(' ')
    name_field = "{}. {}".format(owner[0][0].upper(), ' '.join(owner[1:]))
    if len(name_field) > 18:
        name_field = "{}*".format(name_field[0:17])

    nlabel = input('Enter the number of labels to print: ')
    # Verify that input is of correct type
    try:
        nlabel = int(nlabel)
    except ValueError:
        print("Error: number of labels must be an integer value", \
              file=sys.stderr)
        leave = input('Press any key to exit')
        sys.exit(1)

    description = input("Optional description to be added to the label (enter "
                        "to skip):")
    # Verify that input length is less than 36 characters
    if len(description) > 36:
        print("Error: optional description must be less than 36 characters "
              "long", file=sys.stderr)
        leave = input('Press any key to exit')
        sys.exit(1)

    # Find the sample ID range
    first_id = get_new_id(args.master, args.code, sep=args.sep)
    last_id = first_id + nlabel
    sample_range = "{}-{}".format(first_id, last_id - 1)

    # Write or append to output files
    output = "{}{}_{}_{}.csv".format(args.out_dir, assignee.replace(" ", "_").lower(), \
             sample_range, current_date)

    zpl = "{}{}_{}_{}.zpl".format(args.zpl_dir, assignee.replace(" ", "_").lower(), \
          sample_range, current_date)

    print("Outputting table of sample IDs in {}".format(output), \
          file=sys.stdout)
    nfields = len(args.fields)
    nsep = sep * nfields
    with open(output, 'w') as out_h, open(zpl, 'w') as zpl_h, \
        open(args.master, 'a') as master_h:
        out_h.write('SampleID{0}Owner{0}{1}\n'.format(sep, sep.join(args.fields)))
        for sid in list(range(first_id, last_id)):
            sid_full = str(sid).zfill(args.pad)
            out_h.write('{1}{2}{0}{3}{4}\n'.format(sep, args.code, sid_full, assignee, nsep))
            master_h.write('{1}{2}{0}{3}{4}\n'.format(sep, args.code, sid_full, assignee, nsep))
            zpl_h.write(zpl_template.format(args.code, sid_full, name_field, description))

    # Print labels
    print("Generating {} labels\n".format(str(nlabel)), file=sys.stdout)
#    p = subprocess.Popen([args.net, zpl])
#    try:
#        outs, errs = p.communicate(timeout=60)
#    except TimeoutExpired:
#        proc.kill()
#        outs, errs = p.communicate()

#    if errs:
#        print("{}".format(errs), file=sys.stderr)

    # Allow time to note ID number before exiting
    sleep(2)
    leave = input('Press any key to exit')

if __name__ == "__main__":
    main()
    sys.exit(0)
