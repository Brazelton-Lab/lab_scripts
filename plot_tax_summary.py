#! /usr/bin/env python

"""
Create a krona graph from a taxonomy summary file

Copyright:

    plot_tax_summary  Create a krona graph from a taxonomy summary file

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
__date__ = "2015-03-09"

import sys
import argparse
import os

def argument_parser():
    parser = argparse.ArgumentParser(description="Generate a krona chart from "
                                     "a taxonomy summary file")
    parser.add_argument('summary', metavar="TAX_SUM",
                        help="taxonomy summary file")
    parser.add_argument('pref', metavar="PREFIX",
                        type=str,
                        help="prefix for naming of the krona file")
    parser.add_argument('-s', '--split',
                        action="store_true",
                        help="output into separate files by sample")
    return parser

def file_check(infile, mode):
    try:
        fh = open(infile, mode)
        fh.close()
    except IOError as e:
        print(e)
        sys.exit(1)

def parse_sum(summary, split=False):
    taxonomy = {}
    with open(summary, 'rU') as sum_h:
        header = sum_h.readline().strip().replace(' ', '').split('\t')
        index_map = {}
        for index, col_header in enumerate(header):
            index_map[col_header] = index
        try:
            rank_index = index_map["rankID"]
            tax_index = index_map["taxon"]
            total_index = index_map["total"]
            level_index = index_map["taxlevel"]
            daught_index = index_map["daughterlevels"]
        except KeyError:
            print('Error: incorrect header format')
            sys.exit(1)
        def_cat = [rank_index, tax_index, total_index, level_index, daught_index]
        samples = [(header[i], index_map[header[i]]) for i in 
                  range(len(header)) if i not in def_cat]
        
        taxlevels = []
        for line in sum_h:
            row = line.strip().replace(' ', '').split('\t')
            rank_id = row[rank_index]
            taxon = row[tax_index]
            total = row[total_index]
            taxlevel = row[level_index]
            if taxlevel not in taxlevels:
                taxlevels.append(int(taxlevel))
            taxonomy[rank_id] = {"taxon": taxon, "total": total, "taxlevel": taxlevel}
            if split:
                for sample in samples:
                    sample_name = sample[0]
                    sample_position = sample[1]
                    taxonomy[rank_id][sample_name] = row[sample_position]
    if not split:
        return taxonomy, sorted(taxlevels)[-1]
    else:
        return taxonomy, sorted(taxlevels)[-1], samples

def write_file(outfile, taxonomy, key, max_l):
    with open(outfile, 'w') as out:
        for ident in sorted(taxonomy):
            if int(taxonomy[ident]["taxlevel"]) < max_l or int(taxonomy[ident][key]) == 0:
                continue
            phylogeny = []
            levels = ident.split('.')
            for index in range(len(levels)):
                branch = '.'.join(levels[0: index + 1])
                phylogeny.append(taxonomy[branch]["taxon"]) 
            output = "{} {}\n".format(str(taxonomy[ident][key]), '\t'.join(phylogeny))
            out.write(output)

def main():
    args = argument_parser().parse_args()
    file_check(args.summary, 'rU')
    if args.split:
        taxonomy, max_level, samples = parse_sum(args.summary, args.split)
        for sample in samples:
            abund_key = sample[0]
            outfile = "{}.krona".format(abund_key)
            file_check(outfile, 'w')
            write_file(outfile, taxonomy, abund_key, max_level)
    else:
        taxonomy, max_level = parse_sum(args.summary, args.split)
    
    outfile = args.pref + '.krona'
    file_check(outfile, 'w')
    abund_key = "total"
    write_file(outfile, taxonomy, abund_key, max_level)

if __name__ == "__main__":
    main()
