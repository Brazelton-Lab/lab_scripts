#! /usr/bin/env python
"""
Generates a fasta file of protein sequences and an id mapping file for UniProt
using the m5nr seq and function files.

Usage:

    m5nr_uniprot.py -f out_fasta -m out_mapper sequence_file function_file

"""

from __future__ import print_function

import sys
import os
import argparse
import re
import io
import gzip
import bz2

def open_output(filename):
    """
    Decide how to open an output file for writing based on the file extension
    """
    extension = filename.split('.')[-1]
    if extension == 'gz':
        return gzip.GzipFile(filename, 'wb')
    elif extension == 'bz2':
        return bz2.BZ2File(filename, 'wb')
    else:
        return open(filename, 'w')

def open_input(filename):
    """
    Make a best-effort guess as to how to open the given input file.
    Support for gzip and bzip2 compression.
    """
    file_signatures = {
        "\x1f\x8b\x08": "gz",
        "\x42\x5a\x68": "bz2",
        # "\x50\x4b\x03\x04": "zip"
    }  # Inspired by http://stackoverflow.com/a/13044946/1585509
    try:
        bufferedfile = io.open(file=filename, mode='rb', buffering=8192)
    except IOError:
        print_message("error: unable to open {} for reading".format(filename),
            sys.stderr)
    num_bytes_to_peek = max(len(x) for x in file_signatures)
    file_start = bufferedfile.peek(num_bytes_to_peek)
    compression = None
    for signature, ftype in file_signatures.items():
        if file_start.startswith(signature):
            compression = ftype
            break

    if compression is 'bz2':
        return bz2.BZ2File(filename)
    elif compression is 'gz':
        if not bufferedfile.seekable():
            print_message("gziped data not streamable, pipe through zcat \
                             first", sys.stderr)
        return gzip.GzipFile(filename)
    else:
        return bufferedfile

def main():
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('sequence_file', dest='seq_h',
        type=open_input,
        help="m5nr sequence file [required]")

    parser.add_argument('function_file', dest='func_h',
        type=open_input,
        help="m5nr function file [required]")

    parser. add_argument('-m', '--output-mapper', metavar='FILE',
        dest='mapper_h', type=open_output,
        help="output mapper file name. Support for gzip or bzip2 compression if "
            "the appropriate extension is added to the file name (.gz or .bz2)")

    parser.add_argument('-f', '--output-fasta', metavar='FILE', dest='fasta_h',
        type=open_output,
        help="output fasta file name. Support for gzip or bzip2 compression if "
            "the appropriate extension is added to the file name (.gz or .bz2)")

    args = parser.parse_args()

    records = {}

    r = re.compile("((?<=\()EC.+(?=\)))")
    for line in args.func_h:
        try:
            m5_id, gene_id, product, organism, db = line.strip().split('\t')
        except ValueError:
            raise ""
        matched = r.search(product)
        if not matched:
            ec = ''
        else:
            ec = matched.group()
            product = product.replace('({})'.format(ec), '').rstrip()
        if m5_id not in records:
            records[m5_id] = {'function': [(gene_id, product, ec, organism)]}
        else:
            records[m5_id]['function'].append((gene_id, product, ec, organism))

    for line in args.seq_h:
        m5_id, seq = line.strip().split('\t')
        try:
            seqs[m5_id]['sequence'] = seq
        except KeyError:
            print("no match for id: {}".format(m5_id))
            sys.exit(1)

    for gene in records:
        seq = records[gene]['sequence']
        for function_info in records[gene]['function']:
            ident, product, ec, taxon = function_info
            output = ">{} {}~~~{}~~~{}\n{}\n".format(ident, ec, gene, product, seq)
            fasta_h.write(output)

            gene_len = 3 * len(seq)
            output = "{}\t{}\t{!s}\t{}\n".format(gene, ident, gene_len, taxon)
            mapper_h.write(output)

if __name__ == "__main__":
    main()
    sys.exit(0)
