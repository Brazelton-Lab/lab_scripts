#! /usr/bin/env python

"""
Extract Phylosift output of taxa_summary.py by bin

Copyright:

    phylosift_bins  Extracts Phylosift output of taxa_summary.py by bin

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
from bio_utils.iterators import fasta_iter
from collections import defaultdict
import os
import pysam
import sys

__author__ = 'Alex Hyer'
__email__ = 'theonehyer@gmail.com'
__license__ = 'GPLv3'
__maintainer__ = 'Alex Hyer'
__status__ = 'Production'
__version__ = '1.`.0'


def main(args):
    """Run main program

    Args:
        args (NameSpace): argparse NameSpace of variables influencing program
    """

    # Read and index taxonomy file for random access
    file_index = defaultdict(list)
    args.taxonomy.readline()
    location = args.taxonomy.tell()
    line = args.taxonomy.readline()
    while line:
        read_name = line.strip().split('\t')[0]
        read_name = read_name.split(' ')[0]
        file_index[read_name].append(location)
        location = args.taxonomy.tell()
        line = args.taxonomy.readline()

    for fasta in args.fasta:
        fasta_path = os.path.abspath(fasta)
        out_file = os.path.basename(fasta_path) + '.taxa_summary.txt'
        out_name = os.path.join(args.output_dir + os.sep + out_file)
        out_name = os.path.normpath(out_name)
        unique = []
        with open(out_name, 'w') as out_handle:
            out_handle.write('#Sequence_ID\tHit_Coordinates\t'
                             'NCBI_Taxon_ID\tTaxon_Rank\tTaxon_Name\t'
                             'Cumulative_Probability_Mass\t'
                             'Markers_Hit{0}'.format(os.linesep))
            with open(fasta_path, 'r') as fasta_handle:
                for entry in fasta_iter(fasta_handle):
                    if entry.id in unique:
                        continue
                    else:
                        unique.append(entry.id)
                    for read in args.bam.fetch(entry.id):
                        try:
                            locations = file_index[read.query_name]
                            for location in locations:
                                args.taxonomy.seek(location)
                                out_handle.write(args.taxonomy.readline())
                        except KeyError:
                            continue


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
