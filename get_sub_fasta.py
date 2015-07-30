#!/usr/bin/env python

"""extracts entries from samples based on IDs and writes them to a new file

Usage:

    get_sub_fasta.py <output> [--fastq] <tool> <ids/input> <ids/input>

Synopsis:

    Scans FASTA or FASTQ files for headers that match at least one
    ID. If a header matches an ID, the entire entry is written to
    a new FASTA or FASTQ file.

Required Arguments:

    output:        Name of file to write matches to
    tool:          Which search alorithm to use
    ids:           One or more IDs to search for (see Tools below)
                   IDs can be written as regular expressions
    input:         One or more FASTA or FASTQ files to search (see Tools below)

Optional Arguments:

    --fastq:       Specifies input files as FASTQ files [Default: FASTA]

Tools:

    multi_id:      This tools searchs a single FASTA or FASTQ file with
                   multiple IDs. If a header matches any of the IDs, it is
                   written to the output file. If specified, the last two
                   required arguments are: <input> <ids>

    multi_sample:  This tool searchs multiple FASTA or FASTQ files with a
                   single ID. All matches are written to the same output file.
                   If specified, the last two required arguments are:
                   <ids> <input>
"""

from __future__ import print_function
import argparse
import re
from screed.fasta import fasta_iter
from screed.fastq import fastq_iter
import sys

__author__ = 'Alex Hyer'
__version__ = '0.0.0.1'


def compile_ids(ids):
    """Converts IDs to raw strings, compiles them, and returns them as a list"""

    compiled_ids = []
    for id in ids:
        compiled_ids.append(re.compile(repr(id)))
    return compiled_ids


def fastaq_iter(file_handle, fastq=False):
    """Yields entries from a FASTA or FASTQ file, recreates original header"""    
 
    fastaq_iterator = fasta_iter if not fastq else fastq_iter
    for entry in fastaq_iterator(file_handle):
        if not fastq:
            entry['name'] = '{0} {1}'.format(entry['name'], \
                                             entry['description]')
        else:
            entry['name'] = '{0} {1}'.format(entry['name'], \
                                             entry['annotations'])
        yield entry
            

def multi_id_extraction(ids, input, fastq=False):            
    """Extracts entries with various IDs from a single sample"""

    entries = []
    compiled_ids = compile_ids(ids)
    with open(input, 'rU') as in_handle:
        for entry in fastaq_iter(in_handle, fastq=fastq):
            for compiled_id in compiled_ids:
                if len(compiled_id.findall(entry['name'])) == 1:
                    to_return = ''
                    if not fastq:
                        to_return = '>{0}\n{1}\n'.format(entry['name'], \
                                                         entry['sequence'])
                    else:
                        to_return = '@{0}\n{1}\n+\n{2}\n'.format(entry['name'],
                                                          entry['sequence'],
                                                          entry['accuracy'])      
                    entries.append(to_return)
                    break
    return entries


def multi_sample_extraction(id, input, fastq=False):
    """Extract entries with a single ID from multiple samples"""

    entries = []
    compiled_id = compile_ids([id])[0]
    for file in input:
        with open(file, 'rU') as in_handle:
            for entry in fastaq_iter(in_handle, fastq=fastq):
                if len(compiled_id.findall(entry['name'])) == 1:
                    to_return = ''
                    if not fastq:
                        to_return = '>{0}\n{1}\n'.format(entry['name'], \
                                                         entry['sequence'])
                    else:
                        to_return = '@{0}\n{1}\n+\n{2}\n'.format(entry['name'],
                                                          entry['sequence'],
                                                          entry['accuracy'])      
                    entries.append(to_return)
    return entries
                    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.
                                     RawDescriptionHelpFormatter)
    parser.add_argument('output',
                        default=None,
                        nargs='?',
                        help='name of FASTA or FASTQ file to write')
    parser.add_argument('--fastq',
                        default=False,
                        help='specify input file as FASTQ file [Default: FASTA]')

    tools = parser.add_subparsers(dest='tool',
                                  help='program options')

    multi_sample = tools.add_parser('multi_sample',
                                    help='extract one ID from ' \
                                         'multiple samples')
    multi_sample.add_argument('ids', metavar='ID',
                              default=None,
                              nargs='?',
                              help='ID to extract from samples')
    multi_sample.add_argument('input',
                              default=None,
                              nargs='*',
                              help='FASTA or FASTQ files to extract ids from')

    multi_id = tools.add_parser('multi_id',
                                help='extract multiple IDS from ' \
                                     'one sample') 
    multi_id.add_argument('input',
                          default=None,
                          nargs='?',
                          help='FASTA or FASTQ file to extract IDs from')
    multi_id.add_argument('ids', metavar='IDs',
                          default=None,
                          nargs='*',
                          help='IDs to extract from the sample')

    args = parser.parse_args()

    if args.output is None:
        print(__doc__)
        sys.exit(0)
    elif args.tool == 'multi_sample':
        entries = multi_sample_extraction(args.ids, args.input, \
                                          fastq=args.fastq)
    elif args.tool == 'multi_id':
        entries = multi_id_extraction(args.ids, args.input,
                                      fastq=args.fastq)
    else:
        print('Must specify multi_sample or multi_id mode.')
        sys.exit(1)

    with open(args.output, 'w') as out_handle:
        for entry in entries:
            out_handle.write(entry)

    sys.exit(0)

