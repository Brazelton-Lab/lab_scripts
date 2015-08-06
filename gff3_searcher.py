#!/usr/bin/env python

"""
"""

import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                    formatter_class=argparse.
                                    RawDescriptionHelpFormatter)
    parser.add_argument('-f', '--fields',
                        nargs='*',
                        choices=[
                            'all',
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
                        help='GFF3 files to search for IDs ' \
                             '[Default: config file]')
    parser.add_argument('--ids', metavar='IDs',
                        nargs='*',
                        help='IDs to search GFF3 files for')
    parser.add_argument('-o', '--output',
                        action='store_true',
                        help='output file to write findings to ' \
                             '[Default: STDOUT]')
    parser.add_argument('-of', '--output_format', metavar='Output Format',
                        choices=[
                            'fasta',
                            'gff3',
                            'tsv'
                        ],
                        default='fasta',
                        help='Output Format to write the Output File in ' \
                             '[Default: fasta]')
    args=parser.parse_args()
