#!/usr/bin/env python3

"""
gff3_searcher v. 1.5.0a1 - a program to filter annotations

Usage:

    gff3_searcher.py [--coverage] [--exact] [--fields] [--gff3_files]
                     [--line_search] [--ids] [--output_dir]
                     [--output_format] [--whole-contig]

Synopsis:

    gff3_searcher searches the various fields of a GFF3 file for a match
    to any of the user-specified IDs. A FASTA or GFF3 file is then written
    where each entry matched at least one ID.

Optional Arguments:

    --coverage      Read a table containing coverage data on the contigs
                    annotated in the GFF3 file. If the output format is FASTA,
                    then specifying this option will add the coverage data
                    to the FASTA headers or GFF3 attributes column.
    --exact         If this flag is specified, --flags must also be specified.
                    If present, hits will only be written of a field perfectly
                    matches an ID.
    --fields        Specifies which fields of the GFF3 file to search
                    for IDs in. By default, all fields are searched.
                    More information on the different fields is found
                    below.
    --gff3_files    GFF3 files to search. By default, a config file is read
                    for a default directory and all GFF3 files in that
                    directory are searched.
    --ids           One or more IDs to search for in the GFF3 files.
                    The IDs can be gene names, keywords, gene identifiers,
                    etc. The IDs are not case-sensitive. Each ID must be
                    seperated by a space. If a multi-word ID is desired,
                    wrap it in quotations, i.e. "ATP synthase" searches for
                    IDs containing the whole term "ATP synthase". If this
                    argument is left out, gff3_searcher will read IDs from
                    stdin line by line.
    --line_search   Search entire line of GFF3 fie instead of attributes
    --ouput_dir     Directory to write output files to. Default is the
                    current working directory. Output file names are as
                    follows: [gff3_file].hits.[fasta/gff]
    --output_format The format to write the output files in. By default,
                    output is in GFF3 format. Currently, the only other
                    option is FASTA format.
    --whole_contig  When FASTA output is specified, this argument retrieves the
                    whole contig rather than just the annotated sequence. Since
                    multiple proteins ca occur per contig, the annotation data
                    is removed from the FASTA header.

Fields:

        PROKKA 1.12-beta (and likely other versions of PROKKA) have a very
    specific format for the last column of a GFF3 file (the attributes column).
    Below is an explanation of the fields that gff3_searcher can specifically
    search (as opposed to searching the whole attributes column)*:

    id:           The unique PROKKA given ID to an annotation. This field is
                  only useful to search when one wants to specifically pull
                  out the information on a PROKKA entry of known ID.
    gene:         The short gene name, i.e. rnfC, minD_1, etc.
    inference:    The protein in the database that PROKKA was run on that the
                  annotation was similar to, aka if PROKKA was run on the
                  UniProt database (default) then this field will contain
                  the UniProt ID that the annotation had a high similarity to.
    locus_tag:    Exact same as id above.
    product:      The 'common name' of the annotation based on the inference.
                  I.e. Heptaprenyl diphosphate synthase component I,
                  Rod shape-determining protein MreB, etc. This field is the
                  one field that is quite general and acts best like a
                  'search engine'. For example, if one of the IDs to search
                  with is 'meth', all (or more safely stated, nearly all)
                  genes related to methane use will be extracted.
    ec_number:    The Enzyme Commision Number, if one is available.

    *Note: gff3_searcher now takes arbitrary field names but these fields are
    often the most useful.

Copyright:

    gff3_searcher  searches the various fields of a GFF3 file for a match
    to any of the user-specified IDs

    Copyright (C) 2016  William Brazelton

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
from bio_utils.iterators import fasta_iter, GFF3Reader
import glob
import os
import re
import sys

__author__ = 'Alex Hyer'
__version__ = '1.5.1'

def compile_ids(ids):
    compiled_ids = []
    for id in ids:
        compiled_ids.append(re.compile(id, flags=re.IGNORECASE))
    return compiled_ids


def gff3_line_by_id_retriever(gff3_handle, ids, fields='all', line=False):
    """ids is list of re.compiled ids"""
    gff_reader = GFF3Reader(gff3_handle)

    for entry in gff_reader.iterate(parse_attr=True):
        if line is True:
            for id in ids:
                if len(re.findall(id, entry.write())) != 0:
                    yield entry
                    break
        elif fields != 'all':
            second_loop_break = False
            for field in fields:
                for id in ids:
                    try:
                        if len(re.findall(id, entry.attributes[field])) != 0:
                            yield entry
                            second_loop_break = True
                            break
                    except KeyError:
                        continue
                    if second_loop_break:
                        break
        else:
            for id in ids:
                attr = entry.write().strip().split('\t')[-1]
                if len(re.findall(id, attr)) != 0:
                    yield entry


def read_config():
    try:
        with open('/usr/local/etc/gff3_searcher', 'rU') as config_handle:
            default_directory = config_handle.read().strip()
            return default_directory
    except IOError:
        print('Could not find config file at /usr/local/etc/gff3_searcher')
        user_input = raw_input(
            'Enter a default directory for the config file: ')
        try:
            with open('/usr/local/etc/gff3_searcher', 'w') as config_handle:
                config_handle.write(user_input + '\n')
            return user_input
        except IOError:
            print('Could not write to /usr/local/etc/gff3_searcher.')
            print(
                'You probably do not have write permissions to this directory')
            print(
                'Contact your system administrator about creating this file,')
            print('or use the --gff3_files option to manually specify files.')


def read_conversion_table(conversion_table_file, database):
    temp_dict = {}
    with open(conversion_table_file, 'rU') as table_handle:
        first_line = table_handle.readline()
        databases = first_line.lower().strip().split('\t')
        database = database.lower()
        if database in databases:
            column_number = databases.index(database)
        else:
            print('{0} not in {1}. Available databases are:'.format(
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
    parser.add_argument('-c', '--coverage',
                        default=None,
                        help='table containing coverage for contigs in GFF3')
    parser.add_argument('--exact',
                        action='store_true',
                        help='only write hits that matche an entire ID')
    parser.add_argument('-f', '--fields',
                        nargs='*',
                        default='all',
                        help='Which fields of the GFF3 files to search'
                             'for IDs [Default: all]')
    parser.add_argument('--gff3_files', metavar='GFF3 Files',
                        nargs='*',
                        default=None,
                        help='GFF3 files to search for IDs '
                             '[Default: config file]')
    parser.add_argument('--ids', metavar='IDs',
                        nargs='*',
                        default=sys.stdin,
                        help='IDs to search GFF3 files for')
    parser.add_argument('-l', '--line_search',
                        action='store_true',
                        help='search entire line instead of just attributes')
    parser.add_argument('--output_dir', metavar='Output Directory',
                        default='',
                        help='Directory to write output to')
    parser.add_argument('--output_format', metavar='Output Format',
                        choices=[
                            'fasta',
                            'gff3'
                        ],
                        default='gff3',
                        help='Output Format to write the Output File in '
                             '[Default: gff3]')
    parser.add_argument('--whole_contig', '-w',
                        action='store_true',
                        help='retrieve entire contig rather than just the '
                             'annotated sequence')
    args = parser.parse_args()

    if not args.gff3_files:
        default_directory = read_config()
        args.gff3_files = glob.glob(default_directory + '*.gff')

    args.ids = [id.strip() for id in args.ids]
    if args.exact and args.fields != 'all':
        args.ids = ['^' + id + '$' for id in args.ids]
    elif args.exact and args.fields == 'all':
        print('Must specify --fields with --exact')
        sys.exit(1)
    compiled_ids = compile_ids(args.ids)

    for file in args.gff3_files:
        if args.coverage is not None:
            database = file.split(os.sep)[-1].rsplit('.', 1)[0]
            coverages = read_conversion_table(args.coverage, database)
        with open(file, 'rU') as gff3_handle:
            hits = []
            for line in gff3_line_by_id_retriever(gff3_handle, compiled_ids,
                                                  fields=args.fields,
                                                  line=args.line_search):
                hits.append(line)

            if args.output_format == 'fasta':
                output = file.replace('.gff', '.hits.fasta')
                output_dir = os.getcwd() + os.sep + args.output_dir
                if args.output_dir is not None:
                    if not os.path.isdir(output_dir):
                        os.mkdir(output_dir)
                        output = output.split(os.sep)[-1]
                        output = output_dir + os.sep + output
                    else:
                        output = output.split(os.sep)[-1]
                        output = os.getcwd() + os.sep + output
                with open(output, 'w') as out_handle:
                    for entry in fasta_iter(gff3_handle):
                        for hit in hits:
                            if entry.id == hit.seqid:
                                start = int(hit.start) - 1
                                end = int(hit.end) - 1
                                hit_sequence = entry.sequence[start:end]
                                try:
                                    hit_name = '{0} start_{1} end_{2} ' \
                                               'strand_\"{3}\" ' \
                                               'annotation_\"{4}\"'\
                                               .format(entry.id,
                                                       str(hit.start),
                                                       str(hit.end),
                                                       hit.strand,
                                                       hit.attributes[
                                                                    'product'])
                                except KeyError:
                                    continue
                                if args.whole_contig:
                                    hit_name += ' whole_contig'
                                    hit_sequence = entry.sequence
                                if args.coverage is not None:
                                    rpk = coverages[hit.seqid]
                                    hit_name = hit_name + ' coverage_' + rpk
                                hit_entry = '>{0}\n{1}\n'.format(hit_name,
                                                                 hit_sequence)
                                out_handle.write(hit_entry)
                                if args.whole_contig:
                                    break

            if args.output_format == 'gff3':
                output = file.replace('.gff', '.hits.gff')
                output_dir = os.getcwd() + os.sep + args.output_dir
                if args.output_dir is not None:
                    if not os.path.isdir(output_dir):
                        os.mkdir(output_dir)
                    output = output.split(os.sep)[-1]
                    output = output_dir + os.sep + output
                else:
                    output = output.split(os.sep)[-1]
                    output = os.getcwd() + os.sep + output
                with open(output, 'w') as out_handle:
                    out_handle.write('##gff-version 3\n')
                    for hit in hits:
                        hit_entry = hit.write()
                        if args.coverage is not None:
                            rpk = coverages[hit.seqid]
                            hit_entry = hit_entry.replace('\n', \
                                                          ';coverage:{0}\n'.format(
                                                              rpk))
                        out_handle.write(hit_entry)

    sys.exit(0)
