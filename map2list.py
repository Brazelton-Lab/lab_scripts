#! /usr/bin/env python

"""

Copyright:

    map2list  converts a OTU map file to a list

    Copyright (C) 2016  William Brazelton, Christopher Thornton

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

__author__ = "Christopher Thornton"
__date__ = "2014-12-18"

import sys
import os
import argparse
import re

def argument_parser():
    parser = argparse.ArgumentParser(description="Convert an OTU map file \
             from QIIME to a Mothur-formatted list file")
    parser.add_argument('otu_map',
                        type=str,
                        help="Input OTU map file (the output from QIIME's \
                        pick_otus.py script)")
    parser.add_argument('-c', '--count',
                         type=str,
                         dest="count_file",
                         help="Optional Mothur-formatted count file. This will \
                         produce a new count file containing only the \
                         sequences that are also found in the output list file \
                         (useful for further processing in Mothur)")
    parser.add_argument('-i', '--identity',
                        default=0.97,
                        help="The sequence similarity threshold used in \
                        clustering, a number or none [default 0.97]")
    return parser

def parse_map_file(infile):
    mapper = {}
    with open(infile, 'rU') as in_h:
        for line in in_h:
            array = line.strip().split('\t')
            otu = array[0]
            seqs = array[1:]
            mapper[otu] = seqs
    return mapper

def parse_count_file(infile):
    header = None
    counts = {}
    with open(infile, 'rU') as in_h:
        for line in in_h:
            if "total" in line or "Representative_Sequence" in line:
                header = line
            else:
                seq_id = line.strip().split('\t')[0]
                counts[seq_id] = line
    return header, counts

def sort_by_num(unsorted_list):
    r = re.compile("([a-zA-Z]+)([0-9]+)")
    sorted_tuples = sorted([r.match(i).groups() for i in unsorted_list], \
                           key=return_last)
    sorted_list = [''.join(i) for i in sorted_tuples]
    return sorted_list

def return_last(last):
    return int(last[-1])

def file_check(check_file, mode):
    try:
        fh = open(check_file, mode)
        fh.close()
    except IOError as e:
        print(e)
        sys.exit(1)

def main():
    args = argument_parser().parse_args()
    # verify that input files have valid format
    otu_file = os.path.basename(args.otu_map)
    if otu_file.endswith('.txt'):
        file_check(args.otu_map, 'rU')
        out_list = '.'.join(otu_file.split('.')[:-1]) + '.list'
    else:
        print("{0} does not seem to be a valid OTU map file. Please verify "\
              "{0} is\nformatted correctly and has a valid file extension"\
              .format(otu_file))
        sys.exit(1)
    file_check(out_list, 'w')
    otu_mapper = parse_map_file(args.otu_map)
    if args.count_file:
        file_check(args.count_file, 'rU')
        out_count = '.'.join(otu_file.split('.')[:-1]) + '.count_table'
        file_check(out_count, 'w')
        count_header, seq_counts = parse_count_file(args.count_file)
        with open(out_count, 'w') as count_h:
            count_h.write(count_header)
            for seq_ids in otu_mapper.values():
                for seq_id in seq_ids:
                    if seq_id in seq_counts:
                        count_h.write(seq_counts[seq_id])
    sorted_otus = sort_by_num(otu_mapper.keys())
    header = "label\tnumOtus\t" + '\t'.join(sorted_otus)
    num_otus = len(otu_mapper)
    if type(args.identity) == (type(int()) or type(float())):
        print(type(args.identity))
        ident = '{:.2f}'.format(1 - args.identity)
    else:
        ident = str(args.identity)
    sorted_seqs = [','.join(otu_mapper[i]) for i in sorted_otus]
    with open(out_list, 'w') as list_h:
        list_h.write(header + '\n')
        output = ident + '\t' + str(num_otus) + '\t' + \
                 '\t'.join(sorted_seqs)
        list_h.write(output)
    sys.exit(0)

if __name__ == "__main__":
    main()
