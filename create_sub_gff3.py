#!/usr/bin/env python

"""
extracts entries from samples based on IDs and writes them to a new file

Usage:

    create_sub_fasta.py <output> [--ids] [--files] [--coverage] [--fastq]

Synopsis:

    Scans FASTA or FASTQ files for headers that match at least one
    ID. If a header matches an ID, the entire entry is written to
    a new FASTA or FASTQ file.

Required Arguments:

    output:        Name of file to write matches to
    --ids:         One or more IDs to search for
    --files:       One or more FASTA or FASTQ files to search

Optional Arguments:

    --coverage:    Obtain coverages for the obtained FASTA or FASTQ files.
                   The syntax for this option is as follows:
                   [coverage_table] [GFF file for first input file]
                   [GFF file for second input file] ...
                   WARNING: This is a very specific usage case made for
                   a very specific request. It will NOT give you entry
                   coverage in general. For more info, see Alex Hyer.
    --fastq:       Specifies input files as FASTQ files [Default: FASTA]

Copyright:

    create_sub_gff3.py Create subset GFF file from IDs
    Copyright (C) 2016  William Brazelton, Alex Hyer

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
import argparse
from bio_utils.iterators import GFF3Reader
import re
from screed.fasta import fasta_iter
from screed.fastq import fastq_iter
import sys

__author__ = 'Alex Hyer'
__version__ = '0.0.0.5'


def append_coverages(entries, table, gff3_files):
    return_entries = []
    for file in gff3_files:
        id_converion_dict = create_id_conversion_dict(file)
        database = 'Subbio' + gff3_file.split('.')[0].split(r'/')[-1]
        conversion_table = read_conversion_table(table, database)
        for entry in entries:
            prokka_id = entry.split('\n')[0].split('>')[-1].split(' ')[0]
            if prokka_id in id_conversion_dict:
                contig_id = id_conversion_dict[prokka_id]
                if contig_id in conversion_table:
                    coverage = conversion_table[contig_id]
                    entry_name = entry.split('\n', 1)[0]
                    entry_rest = entry.split('\n', 1)[1]
                    entry_name = entry_name + ' ' + coverage + '\n'
                    new_entry = entry_name + entry_rest
                    return_entries.append(new_entry)
    return return_entries


def create_id_conversion_dict(gff3_file):
    temp_dict = {}
    with open(gff3_file, 'rU') as gff3_handle:
        gff_reader = GFF3Reader(gff3_handle)

        for entry in gff_reader.iterate():
            contig_id = entry.seqid
            prokka_id = entry.attributes['ID']
            if prokka_id not in temp_dict:
                temp_dict[prokka_id] = contig_id

    return temp_dict


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
                                             entry['description'])
        else:
            entry['name'] = '{0} {1}'.format(entry['name'], \
                                             entry['annotations'])
        yield entry


def extract_ids(ids, files, fastq=False):
    """Extract entries with multiple IDs from multiple samples"""

    entries = []
    compiled_ids = compile_ids(ids)
    for file in files:
        with open(file, 'rU') as in_handle:
            for entry in fastaq_iter(in_handle, fastq=fastq):
                for compiled_id in compiled_ids:
                    if len(compiled_id.findall(entry['name'])) == 1:
                        to_return = ''
                        if not fastq:
                            to_return = '>{0}\n{1}\n'.format(entry['name'],
                                                             entry['sequence'])
                        else:
                            to_return = '@{0}\n{1}\n+\n{2}\n'.format(
                                                             entry['name'],
                                                             entry['sequence'],
                                                             entry['accuracy'])
                        # Stop after first ID match
                        break
                        entries.append(to_return)
    return entries


def read_conversion_table(conversion_table_file, database):
    temp_dict = {}
    with open(conversion_table_file, 'rU') as table_handle:
        first_line = table_handle.readline()
        databases = first_line.strip().split('\t')
        if database in databases:
            column_number = databases.index(database)
        else:
            print('{0} not in {1}. Avaialble databases are:'.format(
                  database, conversion_table_file))
            print('\n'.join(databases[1:]))
        for line in table_handle:
            columns = line.strip().split('\t')
            contig_id = columns[0]
            coverage = columns[column_number]
            if coverage != 'None':
                temp_dict[contig_id] = coverage
    return temp_dict


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.
                                     RawDescriptionHelpFormatter)
    parser.add_argument('output',
                        default=None,
                        nargs='?',
                        help='name of FASTA or FASTQ file to write')
    parser.add_argument('--coverage',
                        default=False,
                        help='obtain coverage for FASTA or FASTQ entries,' \
                             ' see help for more details')
    parser.add_argument('--fastq',
                        default=False,
                        action='store_true',
                        help='specify input file as FASTQ file [Default: FASTA]')
    parser.add_argument('--files',
                        default=None,
                        nargs='*',
                        help='Files to extract IDs from')
    parser.add_argument('--ids', metavar='IDs',
                        default=None,
                        nargs='*',
                        help='IDs to extract from input files')
    args = parser.parse_args()

    if args.output is None:
        print(__doc__)
        sys.exit(0)
    elif args.ids is None and args.files is None:
        print('Must specify one or more IDs and one or more files.')
        sys.exit(1)
    else:
        entries = extract_ids(args.ids, args.files, fastq=args.fastq)
        if not args.coverage is None:
            print('WARNING: See help before using the coverage option')
            table = args.coverage[0]
            gff3_files = args.coverage[1:]
            entries = append_coverages(entries, table, gff3_files)

    with open(args.output, 'w') as out_handle:
        for entry in entries:
            out_handle.write(entry)

    sys.exit(0)

