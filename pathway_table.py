#! /usr/bin/env python
""" 
Produce a table containing the metabolic pathway name, coverage, and relative 
abundance in a HumanN-processed sample.
"""

import sys
import os
import argparse
import re

def io_check(infile, mode='rU'):
    try:
        fh = open(infile, mode)
        fh.close()
    except IOError as e:
        sys.exit(str(e))
    return infile

def parse_file(infile, path_table, type):
    file_base = os.path.basename(infile)
    sample_name = file_base.split('.')[0]
    with open(infile, 'rU') as in_h:
        header = in_h.readline()
        if (type not in header.lower()):
            print("error: incorrect header format in {}".format(file_base))
            sys.exit(1)
        for line in in_h:
            line = line.strip().split('\t')
            path_id = line[0]
            value = line[1]
            try:
                path_table[path_id][sample_name][type] = value
            except KeyError:
                try:
                    path_table[path_id][sample_name]= {type: value}
                except KeyError:
                    path_table[path_id] = {sample_name: {type: value}}
    return path_table

def main():
    path_table = {}
    directory = args.directory
    match_cov = re.compile(".+_04a.+-mpt-.+")
    match_abund = re.compile(".+_04b.+-mpt-.+")
    infiles = []
    samples = []
    for root, dirs, files in os.walk(directory):
        for infile in files:
            sample_name = infile.split('.')[0]
            if match_cov.match(infile) or match_abund.match(infile):
                if args.selected:
                    if sample_name in args.selected:
                        infiles.append(os.path.join(root, infile))
                        if sample_name not in samples:
                            samples.append(sample_name)
                else:
                    infiles.append(os.path.join(root, infile))
                    if sample_name not in samples:
                        samples.append(sample_name)
    if not infiles:
        print("error: unable to locate coverage and abundance files in directory \"{}\"".format(directory))
        sys.exit(1)

    if args.prefix:
        prefix = args.prefix
    else:
        prefix = "pathways"
    out_abund = io_check("{}.abundance.csv".format(prefix), 'w')
    out_cov = io_check("{}.coverage.csv".format(prefix), 'w')

    for infile in infiles:
        if match_cov.match(infile):
            type = "coverage"
        elif match_abund.match(infile):
            type = "abundance"
        path_table = parse_file(infile, path_table, type)

    if args.foam:
        path_map = {}
        with open(foam_file, 'rU') as fh:
            header = fh.readline()
            for line in fh:
                col = line.strip().split('\t')
                path_id = col[0]
                pathway = col[1:]
                path_map[path_id] = ';'.join(pathway)

    cov_h = open(out_cov, 'w')
    abund_h = open(out_abund, 'w')

    header = "Pathway\t{}\n".format('\t'.join(sorted(samples)))
    cov_h.write(header)
    abund_h.write(header)
    for path_id in sorted(path_table):
        out_cov = []
        out_abund = []
        try:
            out_path = path_map[path_id]
        except (KeyError, NameError):
            out_path = path_id
        out_cov.append(out_path)
        out_abund.append(out_path)
        for sample in sorted(samples):
            try:
                coverage = path_table[path_id][sample]["coverage"]
            except KeyError:
                coverage = 0
            try:
                abundance = path_table[path_id][sample]["abundance"]
            except KeyError:
                abundance = 0
            out_cov.append(str(coverage))
            out_abund.append(str(abundance))
        cov_h.write("{}\n".format('\t'.join(out_cov)))
        abund_h.write("{}\n".format('\t'.join(out_abund)))

    cov_h.close()
    abund_h.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Produce a table containing "
        "the metabolic pathways found in a sample through HUMAnN, along with "
        "the coverage and abundance of each pathway")
    parser.add_argument('directory', metavar='DIR', 
                        help="directory containing the coverage and abundance files from HUMAnN "
                        "(files containing \"04a\" and \"04b\")")
    parser.add_argument('-s', '--samples', metavar='PREFIX',
                        nargs='*',
                        dest='selected',
                        help="include only samples with these prefixs. Each "
                        "prefix provided should have a match in the coverage "
                        "and abundance files [default: include all]")
    parser.add_argument('-f', '--foam', metavar='FILE',
                        type=io_check,
                        help="file relating foam pathway ids to gene families")
    parser.add_argument('-p', '--prefix', metavar='PREF',
                        help="prefix for the output files")
    args = parser.parse_args()
    foam_file = args.foam
    main()
    sys.exit(0)

