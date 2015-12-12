#! /usr/bin/env python
"""
For dereplicating sequencing reads. Can check for exact duplicates and 
duplicate reads that are prefix of another. Has support for bz2 and gzip 
compression.
"""
from __future__ import print_function
from __future__ import division

__author__ = "Christopher Thornton"
__date__ = "2015-12-07"
__version__ = "0.4.0"

from screed.fastq import fastq_iter
from itertools import izip
import io
import os
import textwrap
import sys
import gzip
import bz2
import argparse
import hashlib

def open_output(filename):
    """Decide how to open an outfile for writing based on file extension"""
    extension = filename.split('.')[-1]
    if extension == 'gz':
        return gzip.GzipFile(filename, 'wb')
    elif extension == 'bz2':
        return bz2.BZ2File(filename, 'wb')
    else:
        return open(filename, 'w')

def open_input(filename):
    """
    Make a best-effort guess as to how to open/parse the given sequence file.

    Deals with .gz, FASTA, and FASTQ records.
    """
    file_signatures = {
        "\x1f\x8b\x08": "gz",
        "\x42\x5a\x68": "bz2",
        # "\x50\x4b\x03\x04": "zip"
    }  # Inspired by http://stackoverflow.com/a/13044946/1585509
    try:
        bufferedfile = io.open(file=filename, mode='rb', buffering=8192)
    except IOError:
        return None
    num_bytes_to_peek = max(len(x) for x in file_signatures)
    file_start = bufferedfile.peek(num_bytes_to_peek)
    compression = None
    for signature, ftype in file_signatures.items():
        if file_start.startswith(signature):
            compression = ftype
            break
    if compression is 'bz2':
        file_handle = bz2file.BZ2File(filename=bufferedfile)
    elif compression is 'gz':
        if not bufferedfile.seekable():
            raise ValueError("gziped data not streamable, pipe through zcat \
                             first")
        file_handle = gzip.GzipFile(filename=filename)
    else:
        file_handle = bufferedfile

    return file_handle
    
def get_iterator(forward_reads, reverse_reads=None):
    """Return object to iterate over"""
    f_iter = fastq_iter(open_input(forward_reads))
    if reverse_reads:
        r_iter = fastq_iter(open_input(reverse_reads))
        return izip(f_iter, r_iter)
    else:
        return f_iter

def replicate_status(query, queried):
    """
    Return the replicate status of a search. A status of zero means not 
    a replicate, one means that the query is a replicate, and two means that 
    the queried is a replicate
    """
    prefix_status = 0
    exact_status = 0
    if query == queried:
        exact_status = 1
    else:
        query_len = len(query)
        queried_len = len(queried)
        if query_len > queried_len:
            if query[:queried_len] == queried:
                prefix_status = 2
        else:
            if query == queried[:query_len]:
                prefix_status = 1
    return (exact_status, prefix_status)

def compare_seqs(query_id, search_id, key, seq_db):
    """Return header of replicate sequence"""
    # check forward read first. If it is a duplicate then check reverse
    query_f, query_r = seq_db[key][query_id]
    search_f, search_r = seq_db[key][search_id]

    fexact, fprefix = replicate_status(query_f, search_f)
    if fexact or fprefix:
        if query_r:
            rexact, rprefix = replicate_status(query_r, search_r)
            if ((fexact == 1) and (rexact == 1)) or \
                ((fexact == 1) and (rprefix == 1)) or \
                ((fprefix == 1) and (rexact == 1)) or \
                ((fprefix == 1) and rprefix == 1):
                return query_id
            elif ((fexact == 1) and (rprefix == 2)) or \
                ((fprefix == 2) and (rexact == 1)) or \
                ((fprefix == 2) and (rprefix == 2)):
                return search_id
            # also implies (fexact/fprefix == 1/2) and (rprefix/rexact== 2/1)
            else:
                return None
        else:
            if (fexact == 1) or (fprefix == 1):
                return query_id
            elif fprefix == 2:
                return search_id

def search_for_duplicates(seq_iter, dups=None, prefix=False, revcomp=False, sub_size=40):
    """Return all records found to be duplicates of another"""
    if not dups:
        dups = []
    md5 = hashlib.md5

    records = {}
    for record in seq_iter:
        paired = True if len(record) == 2 else False
        if paired:
            ident = record[0].name
            fseq, rseq = (record[0].sequence, record[1].sequence)
        else:
            ident = record.name
            fseq, rseq = (record.sequence, '')

        if not prefix:
            key = md5(fseq + rseq).digest()
            if key not in records:
                records[key] = ident
                continue
            else:
                dups.append(ident)
        elif prefix and (fseq < sub_size or rseq < sub_size):
            print(textwrap.fill("warning: some reads are shorter than {!s}bp. "
                "Consider applying a length filter to the dataset"
                .format(sub_size), 79), file=sys.stderr)
            sys.exit(1)
        else:
            key = md5(fseq[:sub_size] + rseq[:sub_size]).digest()
            if key not in records:
                records[key] = {ident: (fseq, rseq)}
                continue
            else:
                search_group = records[key].keys()
                records[key][ident] = (fseq, rseq)

            for search_ident in search_group:
                duplicate = compare_seqs(ident, search_ident, key, records)
                if duplicate:
                    del records[key][duplicate]
                    dups.append(duplicate)
                    break

    return dups

def main():
    parser = argparse.ArgumentParser(description="Remove exact or prefix duplicates from fastq files")
    parser.add_argument('-f', '--forward', dest='in_f', metavar='FASTQ',
        required=True,
        type=str,
        help="paired forward reads in fastq format. Will treat as single-end "
            "reads if the pair is not provided with argument -r/--reverse "
            "[required] ")
    parser.add_argument('-r', '--reverse', dest='in_r', metavar='FASTQ',
        type=str,
        help="paired reverse reads in fastq format")
    parser.add_argument('-o', '--out-forward', dest='out_f', metavar='FASTQ',
        required=True,
        type=open_output,
        help="output forward reads in fastq format [required]")
    parser.add_argument('-l', '--out-reverse', dest='out_r', metavar='FASTQ',
        type=open_output,
        help="output reverse reads in fastq format (use with -r/--reverse)")
    parser.add_argument('-p', '--prefix',
        action='store_true',
        help="replicate can be a 5' prefix of another read")
    parser.add_argument('-m', '--min-size', dest='min_size', metavar='SIZE',
        type=int,
        default=40,
        help="size of the smallest read in the dataset")
    parser.add_argument('-c', '--rev-comp',
        action='store_true',
        help="replicate can be a reverse compliment of another read")
    args = parser.parse_args()

    in_f = args.in_f
    in_r = args.in_r
    out_f = args.out_f
    out_r = args.out_r
    prefix = args.prefix
    substring_size = args.min_size

    if in_r and not out_r:
        print("error: argument -l/--out-reverse required with use of -r/--reverse")
        sys.exit(1)

    print(textwrap.fill("Starting {} with arguments: {}"
        .format(os.path.basename(__file__), ' '.join(sys.argv[1:]), 79)))

    fastq_iterator = get_iterator(in_f, in_r)
    duplicates = search_for_duplicates(fastq_iterator, prefix=prefix, sub_size=substring_size)
    dups_count = len(duplicates)

    print("Number of duplicates: {}".format(dups_count))
    sys.exit(0)
    # write the output
    fastq_iter = get_iterator(in_f, in_r)
    items_count = 0
    for record in fastq_iter:
        items_count += 1
        if in_r:
            header = record[0].name
            seq, rseq = (record[0].sequence, record[1].sequence)
            qual, rqual = (record[0].quality, record[1].quality)
            annot, rannot = (record[0].annotations, record[1].annotations)
        else:
            header = record.name
            seq = record.sequence
            qual = record.quality
            annot = record.annotations

        if header in duplicates:
            continue
        out_f.write("@{}\n{}\n+\n{}\n".format(header, annot, seq, qual))
        if out_r:
            out_r.write("@{} {}\n{}\n+\n{}\n".format(header, rannot, rseq, rqual))

    dups_count = len(duplicates)
    ratio_dups = dups_count / items_count
    print("\nDereplication Complete\n\nReads/Pairs processed:\t{!s}\nDuplicates "
        "found:\t{!s} ({:.2%})\n\n".format(items_count, dups_count, ratio_dups))

if __name__ == "__main__":
    main()
