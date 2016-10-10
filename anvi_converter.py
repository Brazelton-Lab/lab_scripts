#! /usr/bin/env python

"""Convert various file types to TSVs for use with anvi'o

This program has essentially two functions: splitting PROKKA GFF3 files
into a gene_location and genes file to tell anvi'o both where your gene calls
are and what they consist of, and to take FASTA files and produce a tsv for
anvi'o describing which contig is in which bin.
"""

import argparse
from bio_utils.iterators import fasta_iter, gff3_iter
import os
import sys

__author__ = 'Alex Hyer'
__email__ = 'theonehyer@gmail.com'
__license__ = 'GPLv3'
__maintainer__ = 'Alex Hyer'
__status__ = 'Production'
__version__ = '1.1.4b1'


def main(args):
    """Run program

    Args:
         args (NameSpace): ArgParse arguments controlling program flow
    """

    if args.tool == 'bins':
        for fasta in args.FASTA:
            with open(fasta, 'r') as file_handle:
                for entry in fasta_iter(file_handle):
                    file_name = os.path.basename(file_handle.name)
                    args.output.write('{0}\t{1}{2}'.format(entry.id,
                                                           file_name,
                                                           os.linesep))

    if args.tool == 'prokka':
        with open(args.prefix + '.gene_locations.tsv', 'w') as lh, \
                open(args.prefix + '.genes.tsv', 'w') as gh:

            caller_id = 1
            lh.write('gene_callers_id\tcontig\tstart\tstop\t'
                     'direction\tpartial\tsource\tversion{0}'
                     .format(os.linesep))
            gh.write('gene_callers_id\tsource\taccession\t'
                     'function\te_value{0}'.format(os.linesep))
            for entry in gff3_iter(args.GFF3):
                if entry.type == 'CDS' and \
                        'gene' in entry.attributes.keys():

                    print(entry)

                    # Reformat data for gene locations file
                    direction = 'f' if entry.strand == '+' else 'r'
                    program, version = entry.source.split(':')
                    lh.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\tv{7}{8}'
                             .format(str(caller_id),
                                     entry.seqid,
                                     str(entry.start + 1),
                                     str(entry.end - 1),
                                     direction,
                                     '0', program, version, os.linesep))

                    # Reformat data for genes file
                    gh.write('{0}\t{1}\t{2}\t{3}\t{4}{5}'
                             .format(str(caller_id),
                                     entry.source,
                                     entry.attributes['gene'],
                                     entry.attributes['product'],
                                     '0', os.linesep))

                    caller_id += 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.
                                     RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(title='Tool',
                                       dest='tool')

    bins = subparsers.add_parser('bins',
                                 help='read several FASTA files and produce '
                                      'TSV relating entries to file name')
    bins.add_argument('FASTA',
                      nargs='+',
                      type=str,
                      help='list of space-separated FASTA files where each '
                           'file is a bin')
    bins.add_argument('output',
                      type=argparse.FileType('w'),
                      help='output file')

    prokka = subparsers.add_parser('prokka',
                                   help='convert GFF3 file from Christopher '
                                        'Thornton\'s modified version of '
                                        'PROKKA into two TSVs containing '
                                        'called gene locations and their '
                                        'annotations')
    prokka.add_argument('GFF3',
                        type=argparse.FileType('r'),
                        help='GFF3 file to convert')
    prokka.add_argument('prefix',
                        type=str,
                        help='prefix for output files')
    args = parser.parse_args()

    main(args)

    sys.exit(0)
