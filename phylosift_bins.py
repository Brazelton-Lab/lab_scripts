#! /usr/bin/env python

from __future__ import print_function

"""


"""

import argparse
from bio_utils.iterators import fasta_iter
from collections import defaultdict
import os
import pysam
import sys

__author__ = 'Alex Hyer'
__email__ = 'theonehyer@gmail.com'
__license__ = 'GPLv3'
__maintainer__ = 'Alex Hyer'
__status__ = 'Alpha'
__version__ = '0.0.1a1'


def main(args):
    """Run main program

    Args:
        args (NameSpace): argparse NameSpace of variables influencing program
    """

    file_index = defaultdict(int)

    for fasta in args.fasta:
        for entry in fasta_iter(open(fasta, 'r')):
            out_file = os.sep.join(fasta.name.split('.')[:-1])
            out_name = os.path.abspath(args.output_dir + out_file)
            with open(out_name) as out_handle:
                for read in args.bam.fetch(entry.name):
                    pass


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.
                                     RawDescriptionHelpFormatter)
    parser.add_argument('-b', '--bam',
                        type=pysam.AlignmentFile,
                        help='BAM file mapping all reads to bins')
    parser.add_argument('-d', '--output_dir',
                        type=str,
                        help='relative output directory for phylosift files')
    parser.add_argument('-f', '--fasta',
                        nargs='+',
                        type=str,
                        help='list of space-separated FASTA files where each '
                              'file is a bin')
    parser.add_argument('-t', '--taxonomy',
                        type=argparse.FileType('r'),
                        help='sequence_taxa_summary file from phylosift for'
                             'all reads possibly mapping to bins')
    args = parser.parse_args()

    main(args)

    sys.exit(0)
