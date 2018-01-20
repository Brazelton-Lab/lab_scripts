#! /usr/bin/env python

"""Creates a subset BAM file from a FASTA file

Notes:
    1. This program will spawn three processes and works most quickly when
       it has access to at least three CPUs/threads.
    2. This program requires you have samtools installed and globally
       accessible.

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
from subprocess import PIPE, Popen
import sys

__author__ = 'Alex Hyer'
__email__ = 'theonehyer@gmail.com'
__license__ = 'GPLv3'
__maintainer__ = 'Alex Hyer'
__status__ = 'Production/Stable'
__version__ = '1.0.0'


def main(args):
    """Run program

    Args:
         args (NameSpace): ArgParse arguments controlling program flow
    """

    # Collect contig IDs
    contigs = {entry.id: '' for entry in fasta_iter(args.fasta)}

    # Start two samtools processes for reading and writing BAM files
    with Popen(['samtools', 'view', '-h', args.bam], stdout=PIPE,
               universal_newlines=True) as sam_in, \
         Popen(['samtools', 'view', '-h', '-b', '-o', args.output],
               stdin=PIPE, universal_newlines=True) as bam_out:

        for line in sam_iter(sam_in.stdout, headers=True):
            
            if type(line) is str:
                line = line.strip()

                # Only write out headers corresponding to FASTA entries
                if line.startswith('@SQ') is True:
                    try:
                        assert contigs[line.split('\t')[1].split(':')[1]] == ''
                        bam_out.stdin.write(line + os.linesep)
                    except KeyError:
                        pass

                # Write all other headers
                elif line.startswith('@') is True:
                    bam_out.stdin.write(line + os.linesep)

            # Write out BAM entries corresponding to FASTA entries
            else:
                try:
                    assert contigs[line.rname] == ''
                    bam_out.stdin.write(line.write())
                except KeyError:
                    pass

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
