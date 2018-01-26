#! /usr/bin/env python

"""
This script modifies the Mothur taxonomy-summary file by including the full 
phylogeny in the taxon field. It also outputs additional files at different
taxonomic levels

Copyright:

    modify_tax_summary  modifies the Mothur taxonomy-summary file to including the full phylogeny in the taxon field.

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
__date__ = "2014-12-10"
__version__ = "0.1"

import os
import sys
import argparse
import textwrap

def split_args(arguments):
    levels = []
    arguments = [i.lstrip() for i in arguments.split(',')]
    for argument in arguments:
        try:
            rank = int(argument)
        except ValueError:
            output = textwrap.fill("Error: '{}' must be an integer ranging "
                  "from 0 to 8".format(argument, script), 79)
            print(output)
            sys.exit(1)
        else:
            if rank not in list(range(0, 9)):
                output = textwrap.fill("Error: {} is not a valid taxonomy "
                    "level".format(argument, script), 79)
                print(output)
                sys.exit(1)
    return arguments

def file_check(infile, mode='rU'):
    try:
        fh = open(infile, mode)
    except IOError as e:
        print(e)
        sys.exit(1)
    else:
        fh.close()
    return infile

def parse_summary(infile):
    taxonomy = {}
    with open(infile, 'rU') as in_h:
        header = in_h.readline().strip().replace(' ', '').split('\t')
        index_map = {}
        for index, col_header in enumerate(header):
            index_map[col_header] = index
        try:
            rank_index = index_map["rankID"]
            taxon_index = index_map["taxon"]
            level_index = index_map["taxlevel"]
        except KeyError:
            output = textwrap.fill("Error: can't parse header", 79)
            print('Error: Please verify that the header is formatted \
                  correctly in the taxonomy summary file')
            sys.exit(1)
        col_order = [level_index, rank_index, taxon_index]
        sample_info = [header[i] for i in range(len(header)) \
                       if i not in col_order]
        new_header = '\t'.join([header[i] for i in col_order] + sample_info)

        for line in in_h:
            line_list = line.strip().replace(' ', '').split('\t')
            rank_id = line_list[rank_index]
            taxon = line_list[taxon_index]
            rank_level = line_list[level_index]
            data = [line_list[i] for i in range(len(line_list)) \
                    if i not in col_order]
            taxonomy[rank_id] = {"taxon": taxon, "level": rank_level, 
                                 "data": data}
    return new_header, taxonomy

def write_output(taxonomy, taxon, out_h):
    output = taxonomy[taxon]["level"] + '\t' + taxon + '\t' + \
        taxonomy[taxon]["phylogeny"] + '\t' + \
        '\t'.join(taxonomy[taxon]["data"]) + '\n'
    out_h.write(output)

def main():
    infile = file_check(args.taxonomy)
    outfile = file_check(os.path.basename(infile) + '.mod', 'w')
    tax_levels = {'1': 'domain', '2': 'phylum', '3': 'class', '4': 'order', 
                  '5': 'family', '6': 'genus', '7': 'species', '8': 'strain'}
    tax_files = []
    for rank in args.rank:
        file_ext = tax_levels[str(rank)]
        tax_file = file_check("{}.{}"
            .format(os.path.basename(infile), file_ext), 'w')
        tax_files.append((tax_file, rank))
    header, taxonomy = parse_summary(infile)

    # edit taxon name to include full taxonomic classifications
    with open(outfile, 'w') as out:
        out.write(header + '\n')
        for taxon in sorted(taxonomy):
            phylogeny = None
            phylogeny_rank = []
            split_rank = taxon.split('.')
            for index in range(len(split_rank)):
                clade = '.'.join(split_rank[0:index + 1])
                phylogeny_rank.append(clade)
            phylogeny = ';'.join([taxonomy[i]["taxon"] \
                for i in phylogeny_rank if i != '0'])
            if phylogeny:
                taxonomy[taxon]["phylogeny"] = phylogeny
            else:
                taxonomy[taxon]["phylogeny"] = taxonomy[taxon]["taxon"]
            write_output(taxonomy, taxon, out)

    # write outfiles for given (or default) levels
    for tax_file in tax_files:
        out_name = tax_file[0]
        rank = int(tax_file[1])
        tax_subset = [i for i in taxonomy \
                      if (len(i.split('.')) - 1) == rank]
        with open(tax_file[0], 'w') as out:
            out.write(header + '\n')
            for taxon in sorted(tax_subset):
                write_output(taxonomy, taxon, out)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Modify Mothur taxonomy-summary file to add parent groups "
            "to taxa names. Can also optionally output files containing only "
            "taxa of levels specified.")
    parser.add_argument('taxonomy', metavar='tax_summary',
        type=file_check,
        help="Mothur-formatted taxonomy summary file")
    parser.add_argument('-l', '--tax-level', metavar='RANK', dest='rank',
        type=split_args,
        default="5,6",
        help="comma-separated list of ranks for which outfiles should be "
            "produced. 1: domain, 2: phylum, 3: class, 4: order, 5: family, "
            "6: genus, 7: species, 8: strain [default: 5,6]")
    args = parser.parse_args()
    script = os.path.basename(__file__)
    main()

sys.exit(0)
