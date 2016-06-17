#! /usr/bin/env python

"""Adjust BLAST results on contigs by contig coverage"""

from __future__ import print_function

__version__ = '0.0.0.8'
__author__ = 'Alex Hyer'

import argparse
from bio_utils.iterators import gff3_iter
from bio_utils.iterators import b6_iter
from bio_utils.iterators import fasta_iter
import sys


def create_id_conversion_dict(gff3_file):
    temp_dict = {}
    with open(gff3_file, 'rU') as gff3_handle:
        for entry in gff3_iter(gff3_handle, parse_attr=True):
            contig_id = entry.seqid
            try:
                prokka_id = entry.attributes['ID']
                if prokka_id not in temp_dict:
                    temp_dict[prokka_id] = contig_id
            except KeyError:
                continue
    return temp_dict


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


def obtain_unique_rpks(blast_table_file, id_conversion_dict, conversion_table,
                       short_read_length, eValue, alignLen):
    already_hit = {}
    with open(blast_table_file, 'rU') as blast_handle:
        for entry in b6_iter(blast_handle):
            if eValue is None or float(entry.evalue) <= eValue:
                if alignLen is None or int(entry.align_len) >= alignLen:
                    prokka_id = entry.subject
                    if prokka_id in id_conversion_dict:
                        contig_id = id_conversion_dict[prokka_id]
                        if contig_id in conversion_table:
                            if prokka_id not in already_hit:
                                rpk = (float(conversion_table[contig_id])
                                       / float(short_read_length)) * 1000
                                reads = (rpk / 1000) * entry.align_len
                                already_hit[prokka_id] = (rpk, reads)
    return already_hit


def obtain_annotations(fasta_file):
    temp_dict = {}
    with open(fasta_file, 'rU') as fasta_handle:
        for entry in fasta_iter(fasta_handle):
            if entry.id not in temp_dict:
                temp_dict[entry.id] = ' '.join(entry.id.split()[1:])
    return temp_dict


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.
                                     RawDescriptionHelpFormatter)
    parser.add_argument('fasta_file',
                        help='the FASTA file used as the BLAST database')
    parser.add_argument('blast_table', metavar='BLAST table',
                        help='BLAST results table, output format 6')
    parser.add_argument('prokka_gff3_file',
                        metavar='PROKKA generated GFF3 File',
                        help='Annotation file from PROKKA for ID conversion')
    parser.add_argument('contig_coverage_table',
                        metavar='Contig Coverage Table',
                        help='A table containing coverage values for contigs'
                             ' in the BLAST database')
    parser.add_argument('database',
                        help='the sample name')
    parser.add_argument('short_read_length',
                        help='length of short reads for sample')
    parser.add_argument('output',
                        help='the output file to be written')
    parser.add_argument('--max_e_value',
                        type=float,
                        default=None,
                        help='maximum E-value to accept')
    parser.add_argument('--min_align_len',
                        type=int,
                        default=None,
                        help='minimum alignment length to accept')
    args = parser.parse_args()

    print(
        'Warning: Calculating reads and RPK assuming a short read length of {}'
            .format(args.short_read_length))
    id_conversion_dict = create_id_conversion_dict(args.prokka_gff3_file)
    conversion_table = read_conversion_table(args.contig_coverage_table,
                                             args.database)
    unique_rpks = obtain_unique_rpks(args.blast_table, id_conversion_dict,
                                     conversion_table, args.short_read_length,
                                     args.max_e_value, args.min_align_len)
    annotations = obtain_annotations(args.fasta_file)

    with open(args.output, 'w') as out_handle:
        out_handle.write('Annotations\tRPK\tReads\n')
        for prokka_id in unique_rpks:
            out_handle.write('{0}\t{1}\t{2}\n'.format(annotations[prokka_id],
                                                      unique_rpks[prokka_id][
                                                          0],
                                                      unique_rpks[prokka_id][
                                                          1]))

    sys.exit(0)
