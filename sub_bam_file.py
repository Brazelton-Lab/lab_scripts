#! /usr/bin/env python

"""Creates a subset BAM file from a FASTA file

Copyright:
    sub_bam_file.py  create a subset BAM file from a FASTA file
    Copyright (C) 2017  William Brazelton, Alex Hyer

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

from __future__ import unicode_literals

import argparse
from bio_utils.iterators import fasta_iter, sam_iter
import os
import pysam
import sys

__author__ = 'Alex Hyer'
__email__ = 'theonehyer@gmail.com'
__license__ = 'GPLv3'
__maintainer__ = 'Alex Hyer'
__status__ = 'Planning'
__version__ = '0.1.0a1'


def main(args):
    """Run program

    Args:
         args (NameSpace): ArgParse arguments controlling program flow
    """

    contigs = [entry.id for entry in fasta_iter(args.fasta)]

    os.system('samtools view -h -o temp.sam {0}'.format(args.bam))

    with open('temp.sam', 'r') as sam_in, open('temp.new.sam', 'w') as sam_out:

        for line in sam_in:

            line = line.strip()

            if line.startswith('@HD') is True:
                sam_out.write(line + os.linesep)

            elif line.startswith('@SQ') is True:
                parts = line.split('\t')
                contig = parts[1].split(':')[1]
                if contig in contigs:
                    sam_out.write(line + os.linesep)

            elif line.startswith('@') is True:
                sam_out.write(line + os.linesep)

            else:
                if line.split('\t')[2] in contigs:
                    sam_out.write(line + os.linesep)

    os.system('samtools view -h -b -o {0} temp.new.sam'.format(args.output))
    os.remove('temp.sam')
    os.remove('temp.new.sam')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.
                                     RawDescriptionHelpFormatter)

    parser.add_argument('fasta',
                        type=argparse.FileType('r'),
                        help='FASTA file to get IDs from')
    parser.add_argument('bam',
                        type=str,
                        help='BAM file to filter')
    parser.add_argument('output',
                        type=str,
                        help='filtered BAM file to output')
    args = parser.parse_args()

    main(args)

    sys.exit(0)

