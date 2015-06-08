#! /usr/bin/env python
"""
This program takes one or more html krona graphs from phylosift and converts 
the information into a table.
"""

from __future__ import print_function
import os
import sys
import argparse
import re
import textwrap

def add_value(table, ident, name, magnitude, order):
    try:
        row_size = len(table[ident]['values'])
    except KeyError as e:
        table[ident] = {'name': name, 'values': []}
        # add zeros to previous samples that did not contain this taxon
        for number in range(1, order):
            table[ident]['values'].append(str(0))
        table[ident]['values'].append(str(magnitude))
    else:
        # don't know why this would happen, but check anyway
        if name != table[ident]['name']:
            print_text("error: different ids used for same taxon")
            sys.exit(1)
        # add sample's abundance value for this taxon
        table[ident]['values'].append(str(magnitude))
    return table

def parse_html(infile, table, order):
    r = re.compile("(?<=name=\")(?P<name>.*?)(?=\").*?(?<=id=)(?P<id>\d+)(?=\").*?(?<=<val>)(?P<value>\d+\.?\d+)(?=</val>).*?>")
    # parse through file until tax information encountered
    with open(infile, 'rU') as in_h:
        for line in in_h:
            if line.startswith('<br>'):
                line = line.strip()
                break
    # remove non-informative parts of line
    start_index = line.index('<node')
    nodes = line[start_index:]
    # break into segments for recursion
    nodes = nodes.split('<node ')[1:]

    full_ident = [] # to keep track of parent tax ids
    used = [] # to keep track of the taxa encountered in this iteration
    max_level = 0
    for node in nodes:
        match = r.search(node)
        if not match:
            print("error: format incorrect for {}".format(os.path.basename(infile)))
            sys.exit(1)
        taxon_name = match.group('name')
        ncbi_id = match.group('id')
        full_ident.append(ncbi_id)
        if len(full_ident) > max_level:
            max_level = len(full_ident)
        taxon_id = '.'.join(full_ident)
        used.append(taxon_id)
        taxon_mag = match.group('value')
        table = add_value(table, taxon_id, taxon_name, taxon_mag, order)
        # determine how far up the tax tree to ascend
        end_tags = node.count('</node>')
        if end_tags:
          full_ident = full_ident[: -(end_tags)]

    print(str(max_level - 1))

    # at each level, create unassigned group if the sum of immediate 
    # subroups do not match the magnitude of the group
    unassigned = []
    for parent_id in used:
        subs = []
        parent_level = len(parent_id.split('.'))
        parent_mag = float(table[parent_id]['values'][order - 1])
        sub_ids = [i for i in used if ((len(i.split('.')) == parent_level + 1) and (parent_id.split('.') == i.split('.')[:parent_level]))]
        for sub_id in sub_ids:
            sub_mag = float(table[sub_id]['values'][order - 1])
            subs.append(sub_mag)
        if not subs:
            continue
        subs_mag = sum(subs)
        if parent_mag > subs_mag:
            unassigned_id = parent_id + '.0'
            unassigned_name = "unassigned " + table[parent_id]['name']
            unassigned_mag = parent_mag - subs_mag
            unassigned.append((unassigned_id, unassigned_mag))
            table = add_value(table, unassigned_id, unassigned_name, unassigned_mag, order)
        else:
           continue
    # add a zero value to unused taxon (because each sample has different 
    # taxonomic profile)
    for taxon_id in table:
        if (taxon_id not in used) and (taxon_id not in [i[0] for i in unassigned]):
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
        file_name = os.path.basename(infile)
        try:
            fh = open(infile)
            fh.close()
        except IOError as e:
            print_text("warning: {}. The content of {} will not be included in "
                      "the output".format(e, file_name))
            continue
        sample = file_name.split('.')[0]
        # for writing header in order to get values to correspond with correct file
        samples.append(sample)
        table = parse_html(infile, table, iteration)
        iteration += 1

    try:
        out_h = open(outfile, 'w')
    except NameError:
        out_h = None

    header = "Level\tID\tName\t{}".format('\t'.join(samples))
    final_output(out_h, header)
    for ident in sorted(table):
        level = len(ident.split('.')) - 1
        if args.level:
            if level != args.level:
                continue
        name = table[ident]['name']
        values = table[ident]['values']
        if args.no_id:
            final_output(out_h, "{}\t{}\t{}".format(str(level), name, '\t'.join([str(i) for i in values])))
        else:
            final_output(out_h, "{}\t{}\t{}\t{}".format(str(level), ident, name, '\t'.join([str(i) for i in values])))

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
                        help="write to outfile [default: write to STDOUT]")
    parser.add_argument('-f', '--force',
                        action='store_true',
                        help="force overwrite of previous existing output file")
    parser.add_argument('--no-id',
                        dest='no_id',
                        action='store_true',
                        help="do not output taxon ids")
    args = parser.parse_args()
    main()
