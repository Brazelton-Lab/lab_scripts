#! /usr/bin/python3
"""Combine multiple assemblies into a single FASTA file. Sequence headers 
are also modified during merging to prevent duplicate headers.

Required input is one or more sequence assemblies in FASTA format.

Copyright:

    combine_assemblies merge multiple assemblies into one
    Copyright (C) 2016  William Brazelton

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
import os
import re
from seq_annot.argparse import *
import sys
import textwrap
from time import time

__author__ = 'Christopher Thornton'
__license__ = 'GPLv3'
__maintainer__ = 'Christopher Thornton'
__status__ = "Alpha"
__version__ = '0.0.2'

def calculate_stats(len_dist):
    """Calculate assembly statistics, including number sequences, N50, length 
    of the largest sequence, and total size in bp
    """
    len_dist_sort = sorted(len_dist, reverse=True)

    n = len(len_dist)
    largest = len_dist_sort[0]
    shortest = len_dist_sort[-1] 
    total_size = sum(len_dist)

    n50 = 0
    half = total_size / 2
    genome = 0
    for size in len_dist_sort:
        genome += size
        if genome > half:
            n50 = size
            break

    return((n, n50, largest, shortest, total_size))

def split_by_gap(sequence, gap_seq):
    """Split sequence into one or more subsequences based on the presence of 
    gaps of given length
    """
    if not gap_seq:
        return([sequence])

    # Find indices where gap string starts
    split_ind = [m.start() for m in re.finditer(gap_seq, sequence)]

    # Split sequence between occurrences of gap string
    start = 0
    subseqs = []
    for ind in split_ind:
        subseq = sequence[start: ind]
        subseqs.append(subseq)

        start = ind
        for base in sequence[start:]:
            if base == 'N':
                start += 1
            else:
                break

    # Process last subsequence
    subseq = sequence[start:]
    subseqs.append(subseq)

    return([i for i in subseqs if i])

def do_nothing(*args, **kwargs):
    pass

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('infiles',
        metavar='in.fasta',
        nargs='+',
        help="input one or more sequence assemblies in FASTA format")
    parser.add_argument('-o', '--out',
        metavar='out.fa',
        dest='outfile',
        action=Open,
        mode='wb',
        default=sys.stdout,
        help="output combined assembly in FASTA format [default: output to "
            "stdout]")
    parser.add_argument('-l', '--log',
        metavar='out.log',
        dest='log',
        action=Open,
        mode='wb',
        default=sys.stdout,
        help="track header name changes")
    parser.add_argument('-n', '--names',
        metavar='NAME [,NAME,...]',
        action=ParseSeparator,
        sep=',',
        help="comma-separated list of names to be used as a component of the "
            "output sequence headers [default: use input file basenames]. The "
            "order that the names are provided should reflect the order of "
            "the input files")
    parser.add_argument('-g', '--gap',
        type=int,
        metavar='VALUE',
        dest='gaps',
        default=0,
        help="split scaffolds into contigs when greater than specified number "
            "of consecutive Ns encountered [default: 0]")
    parser.add_argument('-m', '--min',
        type=int,
        metavar='VALUE',
        dest='min_len',
        default=100,
        help="discard contigs with final length smaller than threshold "
            "[default: 100]")
    args = parser.parse_args()

    # Output run information
    all_args = sys.argv[1:]
    print("{} {!s}".format('combine_assembly', __version__), file=sys.stderr)
    print(textwrap.fill("Command line parameters: {}"\
          .format(' '.join(all_args)), 79), file=sys.stderr)

    # Track program run-time
    start_time = time()

    # Assign variables based on user input
    infiles = args.infiles
    outfile = args.outfile
    log = args.log.write if args.log else do_nothing

    gaps = args.gaps
    min_len = args.min_len

    gap_seq = 'N' * gaps

    titles = args.names if args.names else [os.path.basename(i) for i \
        in infiles]

    totals = 0
    for position, infile in enumerate(infiles):
        totals += 1

        lengths = []
        seq_num = 0
        title = titles[position]
        try:
            with open(infile, 'r') as fasta_h:
                for record in fasta_iter(fasta_h):

                    seq = record.sequence

                    subseqs = split_by_gap(seq, gap_seq)
                    for subseq in subseqs:
                        seq_num += 1

                        seq_len = len(subseq)

                        # Don't output contigs smaller than min_len
                        if seq_len < min_len:
                            continue

                        lengths.append(seq_len)

                        header = "{}_{!s}".format(title, seq_num)
                        output = ">{}\n{}\n".format(header, subseq).encode('utf-8')
                        outfile.write(output)

                        output_log = "{}\t{}\n".format(header, record.id).encode('utf-8')
                        log(output_log)

        except FileNotFoundError:
            print("error: unable to locate {} for reading"\
                .format(data_path), file=sys.stderr)
            sys.exit(1)
        else:
            ncontig, n50, largest, shortest, totalbp = calculate_stats(lengths)
            output = "{}\nNumber of contigs: {}\nLength of longest contig "\
                "(bp): {}\nLength of shortest contig (bp): {}\nTotal size "\
                "(bp):\nN50: {}".format(title, ncontig, largest, shortest, \
                totalbp, n50)
            print(output, file=sys.stderr)

    outfile.close()

    # Calculate and print program run-time info
    end_time = time()
    total_time = (end_time - start_time) / 60.0
    print("", file=sys.stderr)
    print("It took {:.2e} minutes to merge {!s} assemblies"\
          .format(total_time, totals), file=sys.stderr)
    print("", file=sys.stderr)


if __name__ == "__main__":
    main()
    sys.exit(0)
