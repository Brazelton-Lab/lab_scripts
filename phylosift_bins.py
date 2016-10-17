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
__version__ = '0.0.1a5'


def main(args):
    """Run main program

    Args:
        args (NameSpace): argparse NameSpace of variables influencing program
    """

    # Read and index taxonomy file for random access
    file_index = defaultdict(int)
    args.taxonomy.readline()
    location = args.taxonomy.tell()
    for line in args.taxonomy:
        read_name = line.strip().split('\t')[0]
        file_index[read_name] = location
        location = args.taxonomy.tell()

    for fasta in args.fasta:
        with open(fasta, 'r') as fasta_handle
        for entry in fasta_iter(fasta_handle):
                out_file = os.sep.join(fasta.name.split('.')[:-1])
                out_name = os.path.abspath(args.output_dir + out_file)
                with open(out_name) as out_handle:
                    out_handle.write('#Sequence_ID\tHit_Coordinates\t'
                                     'NCBI_Taxon_ID\tTaxon_Rank\tTaxon_Name\t'
                                     'Cumulative_Probability_Mass\t'
                                     'Markers_Hit{0}'.format(os.linesep))
                    for read in args.bam.fetch(entry.name):
                        try:
                            args.taxonomy.seek(file_index[read.query_name])
                        except KeyError:
                            continue
                        out_handle.write(args.taxonomy.readline())


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
