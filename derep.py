#! /usr/bin/env python

from __future__import print_function
from screed.openscreed import open_reader
from screed.fastq import fastq_iter
from screed.fasta import fasta_iter
import argparse
import gzip

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-f', '--forward',
                        type=open_reader,
                        help="")
    parser.add_argument('-r', '--reverse',
                        type=open_reader,
                        help="")
    parser.add_argument('-g', '--gzip',
                        action='store_true',
                        help="")
    args = parser.parse_args()

    
    


if __name__ == "__main__":
    main()
