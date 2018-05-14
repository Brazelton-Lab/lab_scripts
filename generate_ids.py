#! /usr/bin/env python3
"""
Generates sample table and new labels based on the number of desired labels
and the last sample identifier used.

Copyright:

    generate_ids.py generates project IDs with accompanying labels and tables
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
            print("Error: no column named 'SampleID' found in {\n}"\
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
                      "integers\n", file=sys.stderr)
                leave = input('Press any key to exit')
                sys.exit(1)

            # Store ID numbers
            try:
                sample_number = int(sample[len(prefix):])
            except ValueError:
                print("Error: sample IDs must consist of one or more "
                      "alphabetical characters followed by one or more "
                      "integers\n", file=sys.stderr)
                leave = input('Press any key to exit')
                sys.exit(1)
            samples.append(sample_number)

    if len(samples) <= 0:
        new = 1
    else:
        new = sorted(samples)[-1] + 1
    return(new)


def wrap_text(field):
    """Split on word, if necessary, to prevent final line from exceeding the 
    maximum number of characters allowed per line.
    """
    max_char_per_line = 18
    
    desc = field.split(' ')
    for position, word in enumerate(desc):
        current_line = desc[0:position + 1]
        nchar = len(' '.join(current_line))
        if nchar >= max_char_per_line:
            final_line = desc[position:]
            final_line_len = len(' '.join(final_line))
            if final_line_len > max_char_per_line:
                # Final line too long, must split word between two lines
                first_line = desc[0:position]

                # Split word into two
                word_subset = ''
                for index, character in enumerate(word):
                    word_subset += character
                    if len(' '.join(first_line + [word_subset])) >= max_char_per_line - 1:
                        break

                if len(word_subset) < len(word):
                    first_line.append(word_subset + '-')
                    final_line = [word[index + 1:]] + desc[position + 1:]

                desc = first_line + final_line
            else:
                break

    return(' '.join(desc))


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
    parser.add_argument('--reprint',
        action='store_true',
        help="reprint a range of labels")
    args = parser.parse_args()

    sep = args.sep

    max_char_per_line = 18

    zpl_template = ("^XA\n^FO30,35^ATN^FDLost City 2018^FS\n"
                    "^FO280,35^FDMM,A{0}{1}^BQN,2,6,H^FS\n"
                    "^FO265,195^ATN^FB174,2,,C^FD{0}{1}^FS\n"
                    "^FO30,95^ARN^FD{2}^FS\n"
                    "^FO30,135^ARN^FD{0}{1}^FS\n"
                    "^FO30,175^ARN^FB250,2,,L^FD{3}^FS\n^XZ\n")

    # Current date in ISO format
    current_date = date.today().isoformat().replace('-', '')

    # Verify that length of code plus padding is less than max sampleID length
    if len(args.code) + args.pad > 8:
        print("Error: the sample ID contains too many characters for this "
              "label. Please contact an administrator for assistance.")
        leave = input('Press any key to exit')
        sys.exit(1)

    # Verify that network printer batch file exists
    if not os.path.exists(args.net):
        print("Error: cannot find location of the network printer batch file. "
              "Please contact an administrator for assistance")
        leave = input('Press any key to exit')
        sys.exit(1)

    # Make output directories if they don't already exist
    for directory in [args.out_dir, args.zpl_dir]:
        if not os.path.exists(directory):
            os.mkdirs(directory)

    # Check if master table exists and is not empty
    if not os.path.isfile(args.master):
        with open(args.master, 'w') as master_h:
            master_h.write("SampleID{0}Owner{0}{1}\n".format(sep, \
                           sep.join(args.fields)))

    # Obtain input from user
    assignee = input('Enter your name (e.g. Jane Doe): ')
    # Verify that input does not contain separator character
    if args.sep in assignee:
        print("Error: name '{}' cannot contain the same character '{}' that "
              "is used as the field separator\n".format(assignee, sep))
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
        print("Error: number of labels must be an integer value\n", \
              file=sys.stderr)
        leave = input('Press any key to exit')
        sys.exit(1)

    description = input("Enter optional description to be added to the label, "
                        "up to {!s} characters long (press enter to skip): "\
                        .format(max_char_per_line * 2))

    # Verify that input length is less than maximum total characters allowed
    if len(description) >= max_char_per_line * 2:
        print("Error: optional description must be less than {!s} characters "
              "long\n".format(max_char_per_line * 2), file=sys.stderr)
        leave = input('Press any key to exit')
        sys.exit(1)

    # Verify that each line is less than maximim characters allowed per line
    description = wrap_text(description)

    # Find the sample ID range
    if args.reprint:
        first_id = input('Enter the ID number of the first label to print: ')
        try:
            first_id = int(first_id)
        except ValueError:
            print("Error: ID number must be an integer value\n", \
                  file=sys.stderr)
            leave = input('Press any key to exit')
            sys.exit(1)
    else:
        first_id = get_new_id(args.master, args.code, sep=args.sep)

    last_id = first_id + nlabel
    sample_range = "{}-{}".format(first_id, last_id - 1)

    # Write or append to output files
    zpl = "{}{}_{}_{}.zpl".format(args.zpl_dir, assignee.replace(" ", "_").lower(), \
          sample_range, current_date)

    if not args.reprint:
        output = "{}{}_{}_{}.csv".format(args.out_dir, \
                 assignee.replace(" ", "_").lower(), sample_range, current_date)

        print("\nOutputting table of sample IDs in {}".format(output), \
              file=sys.stdout)

        nfields = len(args.fields)
        nsep = sep * nfields
        with open(output, 'w') as out_h, open(zpl, 'w') as zpl_h, \
            open(args.master, 'a') as master_h:
            out_h.write('SampleID{0}Owner{0}{1}\n'.format(sep, \
                        sep.join(args.fields)))

            for sid in list(range(first_id, last_id)):
                sid_full = str(sid).zfill(args.pad)

                out_h.write('{1}{2}{0}{3}{4}\n'.format(sep, args.code, \
                            sid_full, assignee, nsep))
                master_h.write('{1}{2}{0}{3}{4}\n'.format(sep, args.code, \
                               sid_full, assignee, nsep))
                zpl_h.write(zpl_template.format(args.code, sid_full, \
                            name_field, description))
    else:
        with open(zpl, 'w') as zpl_h:
            for sid in list(range(first_id, last_id)):
                sid_full = str(sid).zfill(args.pad)

                zpl_h.write(zpl_template.format(args.code, sid_full, \
                            name_field, description))

    # Print labels
    plural = "s" if nlabel > 1 else ''
    print("\nGenerating {0} label{1}\n".format(str(nlabel), plural), \
          file=sys.stdout)

    p = subprocess.Popen([args.net, zpl])
    try:
        outs, errs = p.communicate(timeout=60)
    except TimeoutExpired:
        proc.kill()
        outs, errs = p.communicate()

    if errs:
        print(errs, file=sys.stderr)

    # Allow time to read output before exiting
    sleep(2)

    if args.reprint:
        print("\nNew labels have been printed. Please edit the sample tables "
              "manually if applicable")

    leave = input('Press any key to exit')

if __name__ == "__main__":
    main()
    sys.exit(0)
