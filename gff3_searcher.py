#!/usr/bin/env python

"""
"""

from __future__ import print_function

import argparse
from bio_utils.iterators.fasta import fasta_iter
from bio_utils.iterators.gff3 import gff3_iter
import glob
import os
import re
import sys


def compile_ids(ids):
    compiled_ids = []
    for id in ids:
        compiled_ids.append(re.compile(id))
    return compiled_ids


def gff3_line_by_id_retriever(gff3_handle, ids, fields='all'):
    """ids is list of re.copmpiled ids"""
    for entry in gff3_iter(gff3_handle, prokka=True):
         if fields != 'all':
             second_loop_break = False
             for field in fields:
                 for id in ids:
                     try:
                         if len(re.findall(id, entry[field])) != 0:
                             yield entry
                             second_loop_break = True
                             break
                     except KeyError:
                         continue
                     if second_loop_break:
                         break
         else:
             for id in ids:
                 if len(re.findall(id, entry['attributes'])) != 0:
                     yield entry
                 break
             


def read_config():
    try:
        with open('/usr/local/etc/gff3_searcher', 'rU') as config_handle:
            default_directory = config_handle.read().strip()
            return default_directory
    except IOError:
        print('Could not find config file at /usr/local/etc/gff3_searcher')
        user_input = raw_input('Enter a default directory for the config file: ')
        try:
            with open('/usr/local/etc/gff3_searcher', 'w') as config_handle:
                config_handle.write(user_input + '\n')
            return user_input
        except IOError:
            print('Could not write to /usr/local/etc/gff3_searcher.')
            print('You probably do not have write permissions to this directory')
            print('Contact your system administrator about creating this file,')
            print('or use the --gff3_files option to manually specify files.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.
                                     RawDescriptionHelpFormatter)
    parser.add_argument('-f', '--fields',
                        nargs='*',
                        choices=[
                            'id',
                            'gene',
                            'inference',
                            'locus_tag',
                            'product',
                            'ec_number'
                        ],
                        default='all',
                        help='Which fields of the GFF3 files to search' \
                             'for IDs [Default: all]')
    parser.add_argument('--gff3_files', metavar='GFF3 Files',
                        nargs='*',
                        default=None,
                        help='GFF3 files to search for IDs ' \
                             '[Default: config file]')
    parser.add_argument('--ids', metavar='IDs',
                        nargs='*',
                        required=True,
                        help='IDs to search GFF3 files for')
    parser.add_argument('--output_directory', metavar='Output Directory',
                        default=None,
                        help='Directory to write output to')
    parser.add_argument('--output_format', metavar='Output Format',
                        choices=[
                            'fasta',
                            'gff3'
                        ],
                        default='fasta',
                        help='Output Format to write the Output File in ' \
                             '[Default: fasta]')
    args=parser.parse_args()

    if not args.gff3_files:
        default_directory = read_config()
        args.gff3_files = glob.glob(default_directory + '*.gff')

    compiled_ids = compile_ids(args.ids)

    for file in args.gff3_files:
        with open(file, 'rU') as gff3_handle:
            hits = []
            for line in gff3_line_by_id_retriever(gff3_handle, compiled_ids, \
                                                  fields=args.fields):
                hits.append(line)

            if args.output_format == 'fasta':
                for entry in fasta_iter(gff3_handle):
                    for hit in hits:
                        if entry['name'] == hit['seqid']:
                            start = int(hit['start']) - 1
                            end = int(hit['end']) - 1
                            hit_sequence = entry['sequence'][start:end]
                            hit_name = '{0} start_{1} end_{2} strand_"{3}" ' \
                                       'annotation_{4}'.format(entry['name'],
                                       hit['start'], hit['end'], hit['strand'],
                                       hit['product'])
                            hit_entry = '>{0}\n{1}\n'.format(hit_name,
                                                             hit_sequence)
                            output = file.replace('.gff', '.hits.fasta')
                            if args.output_directory is not None:
                                output = output.split(os.sep)[-1]
                                output = args.output_directory + output
                            with open(output, 'a') as out_handle:
                                out_handle.write(hit_entry)

            if args.output_format == 'gff3':
                output = file.replace('.gff', '.hits.gff')
                if args.output_directory is not None:
                    output = output.split(os.sep)[-1]
                    output = args.output_directory + output
                with open(output, 'w') as out_handle:
                    out_handle.write('##gff-version 3\n')
                    for hit in hits:
                        hit_entry = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n' \
                                    .format(hit['seqid'], hit['source'],
                                    hit['type'], hit['start'], hit['end'],
                                    hit['score'], hit['strand'], hit['phase'],
                                    hit['attributes'])
                        out_handle.write(hit_entry)

    sys.exit(0)
