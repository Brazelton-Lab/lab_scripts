#!/usr/bin/env python

"""
Separate sequences into different files based on header tags

Usage:
separate_tagged.py -i <interleaved_reads> | 

Copyright:

    separate_seqs  Separate sequences into different files based on header tags

    Copyright (C) 2016  William Brazelton, Christopher Thornton

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.﻿
"""

__version__ = '0.0.1'
__author__ = 'Christopher Thornton'

import argparse
import gzip
import bz2file
from screed.fastq import fastq_iter
from screed.fasta import fasta_iter
from screed.openscreed import open_reader

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    format_group = parser.add_mutually_exclusive_group(required=True)
    format_group.add_argument()
    format_group.add_argument()
    parser.add_argument('-t', '--tag',
                        type=str,
                        help="")
    parser.add_argument('-f', '--filter',
                        type=file_opener,
                        help="")
    args = parser.parse_args()

    

if name == "__main__":
    main()
    sys.exit(0)
