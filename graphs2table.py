#! /usr/bin/env python
"""
Program takes one or more html krona graph files from phylosift and converts
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
            print(name, table[ident]['name'], ident)
            print_text("error: different ids used for same taxon")
            sys.exit(1)
        # add sample's abundance value for this taxon
        table[ident]['values'].append(str(magnitude))
    return table

def parse_html(infile, table, order, max_level=None):
    r = re.compile("(?<=name=\")(?P<name>.*?)(?=\").*?(?<=id=)(?P<id>\d+)(?=\")"
        ".*?(?<=<val>)(?P<value>\d+\.?\d+)(?=</val>).*?>")
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
    for node in nodes:
        match = r.search(node)
        if not match:
            print("error: bad format: {}".format(os.path.basename(infile)))
            sys.exit(1)
        taxon_name = match.group('name')
        ncbi_id = match.group('id')
        full_ident.append(ncbi_id)
        taxon_id = '.'.join(full_ident)
        used.append(taxon_id)
        taxon_mag = match.group('value')
        table = add_value(table, taxon_id, taxon_name, taxon_mag, order)
        # determine how far up the tax tree to ascend
        end_tags = node.count('</node>')
        if end_tags:
          full_ident = full_ident[: -(end_tags)]

    # at each level, create unassigned group if the sum of immediate 
    # subroups do not match the magnitude of the group
    for parent_id in used:
        subs = []
        parent_level = len(parent_id.split('.'))
        parent_mag = float(table[parent_id]['values'][order - 1])
        sub_ids = [i for i in used if ((len(i.split('.')) == parent_level + 1)
            and (parent_id.split('.') == i.split('.')[:parent_level]))]
        for sub_id in sub_ids:
            sub_mag = float(table[sub_id]['values'][order - 1])
            subs.append(sub_mag)
        if not subs:
            continue
        subs_mag = sum(subs)
        if parent_mag > subs_mag:
            unass_id = parent_id + '.0'
            unass_name = "unassigned " + table[parent_id]['name']
            unass_mag = parent_mag - subs_mag
            used.append(unass_id)
            table = add_value(table, unass_id, unass_name, unass_mag, order)

    # create unassigned group if category doesn't terminate at chosen level
    parents = []
    m = re.compile("((^unassigned )?(?P<taxon>.+))")
    for ident in sorted(used, reverse=True, key=sort_by_level):
        level = len(ident.split('.')) - 1
        if max_level != None:
            if level < max_level:
               if ident.split('.') in [i.split('.')[:level + 1] for i in parents]:
                   parents.append(ident)
               else:
                   value = table[ident]["values"][order - 1]
                   name = "unassigned " + m.search(table[ident]['name']).group('taxon')
                   parents.append(ident)
                   ident = "{}.{}".format(ident, '.'.join(["0" for i in range(max_level - level)]))
                   used.append(ident)
                   table = add_value(table, ident, name, value, order)
            elif level == max_level:
                parents.append(ident)

    # add a zero value to unused taxon (because each sample has different 
    # taxonomic profile)
    for taxon_id in table:
        if (taxon_id not in used):
            table[taxon_id]['values'].append("0")
    return table
    
def print_text(text):
    print(textwrap.fill(text, 79))

def final_output(fh, output):
    try:
        fh.write(output + '\n')
    except AttributeError:
        print(output)

def sort_by_level(ident):
    """function to be used as the key for sorting identities"""
    ident = ident.split('.')
    level = len(ident)
    return level

def main():
    user_level = args.level
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
        # for writing header in order to get values to correspond with the
        # correct file
        samples.append(sample)
        table = parse_html(infile, table, iteration, max_level=user_level)
        iteration += 1

    try:
        out_h = open(outfile, 'w')
    except NameError:
        out_h = None

    totals = [0 for i  in range(len(samples))]
    header = "Level\tID\tName\t{}".format('\t'.join(samples))
    final_output(out_h, header)
    for ident in sorted(table):
        level = len(ident.split('.')) - 1
        if level != user_level:
            continue
        name = table[ident]['name']
        values = table[ident]['values']
        final_output(out_h, "{}\t{}\t{}\t{}".format(str(level), ident, name, 
            '\t'.join([str(i) for i in values])))
        for index, abund in enumerate(values):
            if abund == "0":
                continue
            totals[index] += float(abund)
    if args.debug:
        print("sample total(s): {}".format('\t'.join(table["1"]["values"])))
        print("level total(s): {}".format('\t'.join([str(i) for i in totals])))

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
    parser.add_argument('-d', '--debug',
                        action='store_true',
                        help="debug program")
    args = parser.parse_args()
    main()
