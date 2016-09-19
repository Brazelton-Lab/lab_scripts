#! /usr/bin/env python

from __future__ import print_function

"""BLAST sequences from PROKKA annotations


"""

import argparse

__author__ = 'Alex Hyer'
__email__ = 'theonehyer@gmail.com'
__license__ = 'GPLv3'
__maintainer__ = 'Alex Hyer'
__status__ = 'Alpha'
__version__ = '0.0.1a2'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.
                                     RawDescriptionHelpFormatter)
    parser.add_argument('-e', '--e_value',
                        default=10.0,
                        type=float,
                        help='Maximum E-Value of alignment permitted')
    parser.add_argument('-f', '--fna',
                        required=True,
                        type=argparse.FileType,
                        help='FNA file from PROKKA containing nucleotide '
                             'sequences of annotated proteins')
    parser.add_argument('-i', '--id',
                        required=True,
                        type=argparse.FileType,
                        help='line-separated list of IDs designating which '
                             'PROKKA IDs to BLAST with')
    parser.add_argument('-p', '--program',
                        default='blastn',
                        type=str,
                        choices=[
                            'blastn',
                            'blastx',
                            'tblastx'
                        ],
                        help='BLAST+ program to search with')
    parser.add_argument('-t', '--top',
                        default=True,
                        action='store_true',
                        help='only return top BLAST result')

