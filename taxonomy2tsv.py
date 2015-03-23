#! /usr/bin/env python
"""modifies mothur taxonomy file by transfering the last name that is not 
"unclassified" or "uncultured" to "unclassified" or "uncultured" assignment,
also removes numbers in parentheses
"""

from __future__ import print_function
import sys
import argparse
import re
import os

def argument_parser():
    parser = argparse.ArgumentParser(description="Convert a mothur-formatted "
                                     "taxonomy file to tsv. Also removes "
                                     "confidence scores and alters "
                                     "\"unclassified\" and \"uncultured\" "
                                     "assignments to include the last useful "
                                     "assignment in the name")
    parser.add_argument('infile', metavar='TAXONOMY',
                        type=file_check,
                        help="mothur-formatted taxonomy file")
    return parser

def file_check(infile, mode='rU'):
    try:
        fh = open(infile, mode)
        fh.close()
    except IOError as e:
        sys.exit(e)
    return infile

def get_taxonomy(tax, infile):
    r = re.compile("(.*)(\(\d+\))")
    with open(infile, 'rU') as in_h:
        for line in in_h:
            columns = line.strip().split()
            seq_id = columns[0]
            hierarchy = columns[1].split(';')[:-1]
            good_names = []
            for index in range(len(hierarchy)):
                matched = r.match(hierarchy[index])
                if not matched:
                    clade = hierarchy[index]
                else:
                    clade = matched.group(1)
                if (clade == "unclassified" or clade == "uncultured"):
                    prev_good_name = good_names[-1]
                    name = "{}_{}".format(prev_good_name, clade)
                else:
                    name = clade
                    good_names.append(name)
                hierarchy[index] = name
            tax[seq_id] = '\t'.join(hierarchy) 
    return tax

def main():
    args = argument_parser().parse_args()
    infile = args.infile
    outfile = os.path.basename(infile) + '.tsv'
    file_check(outfile, 'w')
    taxonomy = {}
    taxonomy = get_taxonomy(taxonomy, infile)

    with open(outfile, 'w') as out:
        for seq_id in taxonomy:
            output = "{}\t{}\n".format(seq_id, taxonomy[seq_id])
            out.write(output)

if __name__ == "__main__":
    main()
