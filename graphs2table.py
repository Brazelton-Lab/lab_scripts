#! /usr/bin/env python
"""
"""

from __future__ import print_function
import os
import sys
import argparse
import re
import textwrap
from HTMLParser import HTMLParser

def parse_html(infile, table, sample_order):
    r = re.compile("(?<=name=\")(?P<name>.*?)(?=\").*?(?<=id=)(?P<id>\d+)(?=\").*?(?<=<val>)(?P<value>\d+\.?\d+)(?=</val>).*?>")
    # parse through file until tax information encountered
    with open(infile, 'rU') as in_h:
        for line in in_h:
            if line.startswith('<br>'):
                line = line.strip()
                break
    # remove non-useful information from line
    start_index = line.index('<node')
    nodes = line[start_index:]
    # break into segments for recursion
    nodes = nodes.split('<node ')[1:]
    full_ident = [] # to keep track of parent tax ids
    used = [] # to keep track of taxa encountered in this iteration
    for node in nodes:
        match = r.search(node)
        name = match.group('name')
        ncbi_id = match.group('id')
        full_ident.append(ncbi_id)
        ident = '.'.join(full_ident)
        used.append(ident)
        magnitude = match.group('value')
        try:
            row_size = len(table[ident]['values'])
        except KeyError as e:
            table[ident] = {'name': name, 'values': []}
            # add zeros to previous samples that did not contain this taxon
            for number in range(1, sample_order):
                table[ident]['values'].append(str(0))
            table[ident]['values'].append(magnitude)
        else:
            # don't know why this would happen, but check anyway
            if name != table[ident]['name']:
                print_text("error: something weird happened the with names")
                sys.exit(1)
            # add sample's abundance value for this taxon
            table[ident]['values'].append(magnitude)
        # determine how far up the tax tree to ascend
        end_tags = node.count('</node>')
        if end_tags:
          full_ident = full_ident[:-(end_tags)]
    for taxon_id in table:
        if taxon_id not in used:
            table[taxon_id]['values'].append(str(0))
    return table

def print_text(text):
    print(textwrap.fill(text, 79))

def final_output(fh, output):
    try:
        fh.write(output + '\n')
    except AttributeError:
        print(output)

def main():
    if args.out:
        outfile = args.out
        try:
            fh = open(outfile)
        except IOError:
            try:
                fh = open(outfile, 'w')
            except IOError as e:
                print(e)
                sys.exit(1)
            else:
                fh.close()
        else:
            fh.close()
            if not args.force:
                print_text("error: {} already exists in the current working "
                          "directory. Use {} with the --force flag to overwrite"
                          .format(outfile, os.path.basename(__file__)))
                sys.exit(1)

    table = {}
    samples = []
    iteration = 1
    for infile in args.infiles:
        try:
            fh = open(infile)
            fh.close()
        except IOError as e:
            file_name = os.path.basename(infile)
            print_text("warning: {}. The content of {} will not be included in "
                      "the output".format(e, file_name))
            continue
        prefix = file_name.split('.')[0]
        # for writing header in order to get values to correspond with correct file
        samples.append(prefix)
        table = parse_html(infile, table, iteration)
        iteration += 1

    try:
        out_h = open(outfile, 'w')
    except NameError:
        out_h = None

    header = "Level\tID\tName\t{}".format('\t'.join(samples))
    final_output(out_h, header)
    for ident in sorted(table):
        level = len(ident.split('.')) -1
        if args.level:
            if level != args.level:
                continue
        name = table[ident]['name']
        values = table[ident]['values']
        final_output(out_h, "{}\t{}\t{}\t{}".format(str(level), ident, name, '\t'.join(values)))

    try:
        out_h.close()
    except AttributeError:
        sys.exit(0)
    else:
        sys.exit(0)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="create table from kronagraphs")
    parser.add_argument('infiles', metavar='HTML',
                        nargs='+',
                        help="krona html file")
    parser.add_argument('-l', '--level', metavar='VALUE',
                        type=int,
                        help="create table for level \"LEVEL\" [default: all-in-one]")
    parser.add_argument('-o', '--out', metavar='FILE',
                        type=str,
                        help="write to outfile instead of STDOUT")
    parser.add_argument('-f', '--force',
                        action='store_true',
                        help="force overwrite of previous existing output file")
    args = parser.parse_args()
    main()
