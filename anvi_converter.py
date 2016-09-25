#! /usr/bin/env python

from __future__ import print_function

"""Convert various file types to TSVs for use with anvi'o


"""

import argparse
from bio_utils.iterators import fasta_iter, gff3_iter
import sys

__author__ = 'Alex Hyer'
__email__ = 'theonehyer@gmail.com'
__license__ = 'GPLv3'
__maintainer__ = 'Alex Hyer'
__status__ = 'Alpha'
__version__ = '0.0.1a2'


def main(args):
    """Run program

    Args:
         args (NameSpace): ArgParse arguments controlling program flow
    """

    pass


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.
                                     RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(title='Tool',
                                       dest='tool')

    bins = subparsers.add_parser('bins',
                                 help='read several FASTA files and produce '
                                      'TSV relating entries to file name')
    bins.add_argument('FASTA files',
                      dest='fasta_files',
                      type=list,
                      help='list of space-separated FASTA files where each '
                           'file is a bin')
    bins.add_argument('output',
                      type=argparse.FileType('w'),
                      help='output file')

    prokka = subparsers.add_parser('PROKKA',
                                   dest='prokka',
                                   help='convert GFF3 file from Christopher '
                                        'Thornton\'s modified version of '
                                        'PROKKA into two TSVs containing '
                                        'called gene locations and their '
                                        'annotations')
    prokka.add_argument('GFF3',
                        dest='gff3',
                        type=argparse.FileType('r'),
                        help='GFF3 file to convert')
    prokka.add_argument('prefix',
                        type='str',
                        help='prefix for output files')

    sys.exit(0)
