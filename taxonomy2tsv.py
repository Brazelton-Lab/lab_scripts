#! /usr/bin/env python
"""
Modifies a mothur taxonomy file by transfering the last name that is not 
"unclassified" or "uncultured" to "unclassified" or "uncultured" assignment,
also removes numbers in parentheses
"""

from __future__ import print_function
import sys
import argparse
import re
import os

def file_check(infile, mode='rU'):
    try:
        fh = open(infile, mode)
    except IOError as e:
        print(e)
        sys.exit(1)
    else:
        fh.close()
    return infile

def main():
    infile = args.infile
    outfile = file_check(os.path.basename(infile) + '.tsv', 'w')

    r = re.compile("(.*)(\(\d+\))")
    with open(infile, 'rU') as in_h:
        with open(outfile, 'w') as out_h:
            for line in in_h:
                line = line.strip().split('\t')
                seq_id = line[0]
                hierarchy = line[1].rstrip(';').split(';')
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
                output = "{}\t{}\n".format(seq_id, '\t'.join(hierarchy))
                out_h.write(output)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert a mothur-formatted "
        "taxonomy file to tsv. Also removes confidence scores and alters "
        "\"unclassified\" and \"uncultured\" assignments to include the last "
        "useful assignment in the name")
    parser.add_argument('infile', metavar='taxonomy',
        type=file_check,
        help="mothur-formatted taxonomy file")
    args = parser.parse_args()
    main()
