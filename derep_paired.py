#! /usr/bin/env python
"""
For dereplicating paired-end sequencing reads. Can check for exact, 5' prefix, 
and reverse-complement duplicates. Has support for bz2 and gzip compression.
"""
from __future__ import print_function
from __future__ import division

__author__ = "Christopher Thornton"
__date__ = "2015-12-07"
__version__ = "0.7.2"

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
    """Decide how to open an output file for writing based on file extension"""
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

def print_message(message, dest=sys.stdout):
    print(textwrap.fill(message, 79), file=dest)

def get_iterator(forward_reads, reverse_reads=None):
    """Return object to iterate over"""
    f_iter = fastq_iter(open_input(forward_reads))
    if reverse_reads:
        r_iter = fastq_iter(open_input(reverse_reads))
        return izip(f_iter, r_iter)
    else:
        return f_iter

def reverse_complement(forward, reverse):
    compliments = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    forward = ''.join([compliments[i] for i in reverse])
    reverse = ''.join([compliments[i] for i in forward])
    return (forward, reverse)
    
def compare_seqs(query, template):
    """
    Return the replicate status of a search. A status of zero means not 
    a duplicate, one means query is exactly duplicate, two means template is 
    a duplicate, and three means query is a prefix duplicate
    """
    status = 0
    query_len, temp_len= (len(query), len(template))
    if query_len == temp_len:
        if query == template:
            status = 1
    elif query_len > temp_len:
        if query[:temp_len] == template:
            status = 2
    # implies query length < template length
    elif query_len < temp_len:
        if query == template[:query_len]:
            status = 3
    return status

def replicate_status(query_f, query_r, search_f, search_r):
    """Return header of replicate sequence"""
    # check forward read first. If it is a duplicate then check reverse
    fstatus = compare_seqs(query_f, search_f)
    if fstatus:
        rstatus = compare_seqs(query_r, search_r)
        if (fstatus == 1 and rstatus == 1) or \
            (fstatus == 1 and rstatus == 3) or \
            (fstatus == 3 and rstatus == 1) or \
            (fstatus == 3 and rstatus == 3):
            return 'query'
        elif (fstatus == 1 and rstatus == 2) or \
            (fstatus == 2 and rstatus == 1) or \
            (fstatus == 2 and rstatus == 2):
            return 'search'
        else:
            return None
    else:
        return None

def search_for_duplicates(seq_iter, dups=None, prefix=False, revcomp=False, 
    sub_size=8):
    """Return all records found to be duplicates of another"""
    if not dups:
        dups = {}
    md5 = hashlib.md5

    records = {}
    for record in seq_iter:
        ident = record[0].name.split()[0]
        fseq, rseq = (record[0].sequence, record[1].sequence)
        if revcomp:
            fcomp, rcomp = reverse_complement(fseq, rseq)

        if not prefix:
            key = md5(fseq + rseq).digest()
            if key in records:
                dups[ident] = (records[key], 'exact')
                continue
            else:
                if revcomp:
                    compkey = md5(fcomp + rcomp).digest()
                    if compkey in records:
                        dups[ident] = (records[compkey], 'revcomp')
                    else:
                        records[key] = ident
                else:
                    records[key] = ident
        elif prefix and (len(fseq) < sub_size or len(rseq) < sub_size):
            message = ("warning: some reads are shorter than {!s}bp. Consider "
                "applying a length filter to the dataset".format(sub_size))
            print_message(message, dest=sys.stderr)
            sys.exit(1)
        else:
            key = md5(fseq[:sub_size] + rseq[:sub_size]).digest()
            if revcomp:
                compkey = md5(fcomp[:sub_size] + rcomp[:sub_size]).digest()
            if key in records:
                items_to_check = records[key].keys()
                records[key][ident] = (fseq, rseq)
                for search_id in items_to_check:
                    fsearch, rsearch = records[key][search_id]
                    duplicate = replicate_status(fseq, rseq, fsearch, rsearch)
                    if not duplicate:
                        continue
                    else:
                        query_size, search_size = (len(fseq + rseq), 
                            len(fsearch + rsearch))
                        if query_size == search_size:
                            dup_type = 'exact'
                        else:
                            dup_type = 'prefix'

                        if duplicate == 'search':
                            del records[key][search_id]
                            dups[search_id] = (ident, dup_type)
                        else:
                            del records[key][ident]
                            dups[ident] = (search_id, dup_type)
                        break
            elif revcomp and compkey in records:
                items_to_check = records[compkey].keys()
                records[key] = {ident: (fseq, rseq)}
                for search_id in items_to_check:
                    fsearch, rsearch = records[compkey][search_id]
                    duplicate = replicate_status(fcomp, rcomp, fsearch, rsearch)
                    if not duplicate:
                        continue
                    else:
                        query_size, search_size = (len(fcomp + rcomp), 
                            len(fsearch + rsearch))
                        if query_size == search_size:
                            dup_type = 'exact-revcomp'
                        else:
                            dup_type = 'prefix-revcomp'

                        if duplicate == 'search':
                            del records[compkey][search_id]
                            dups[search_id] = (ident, dup_type)
                        else:
                            del records[key][ident]
                            dups[ident] = (search_id, dup_type)
                        break
            else:
                records[key] = {ident: (fseq, rseq)}

    return dups

def main():
    parser = argparse.ArgumentParser(description=__doc__,        
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-f', '--forward', dest='in_f', metavar='FASTQ',
        required=True,
        type=str,
        help="paired forward reads in fastq format [required]")
    parser.add_argument('-r', '--reverse', dest='in_r', metavar='FASTQ',
        required=True,
        type=str,
        help="paired reverse reads in fastq format")
    parser.add_argument('-o', '--out-forward', dest='out_f', metavar='FASTQ',
        required=True,
        type=open_output,
        help="output forward reads in fastq format [required]")
    parser.add_argument('-v', '--out-reverse', dest='out_r', metavar='FASTQ',
        required=True,
        type=open_output,
        help="output reverse reads in fastq format (use with -r/--reverse)")
    parser.add_argument('-l', '--log', metavar='LOG',
        type=open_output,
        help="log file to keep track of duplicates")
    parser.add_argument('-p', '--prefix',
        action='store_true',
        help="replicate can be a 5' prefix of another read")
    parser.add_argument('-m', '--min-size', dest='min_size', metavar='SIZE',
        type=int,
        default=8,
        help="size of the smallest read in the dataset. Should be used with "
            "-p/--prefix [default: 8]")
    parser.add_argument('-c', '--rev-comp', dest='rev_comp',
        action='store_true',
        help="replicate can be the reverse-complement of another read")
    args = parser.parse_args()

    prog_name = os.path.basename(__file__)
    all_args = ' '.join(sys.argv[1:])
    in_f = args.in_f
    in_r = args.in_r
    out_f = args.out_f
    out_r = args.out_r
    log_h = args.log
    prefix = args.prefix
    rev_comp = args.rev_comp
    substring_size = args.min_size

    message = "Starting {} with arguments: {}".format(prog_name, all_args)
    print_message(message)

    seq_iterator = get_iterator(in_f, in_r)
    duplicates = search_for_duplicates(seq_iterator, 
        prefix=prefix, revcomp=rev_comp, sub_size=substring_size)
    dups_count = len(duplicates)

    if log_h:
        log_h.write("Duplicate\tTemplate\tType\n")
        for dup in duplicates:
            template, dup_type = duplicates[dup]
            log_h.write("{}\t{}\t{}\n".format(dup, template, dup_type))

    # write the output
    fastq_iter = get_iterator(in_f, in_r)
    items_count = 0
    for record in fastq_iter:
        items_count += 1
        header, rheader = (record[0].name, record[1].name)
        seq, rseq = (record[0].sequence, record[1].sequence)
        qual, rqual = (record[0].quality, record[1].quality)
        if header.split()[0] in duplicates:
            continue
        out_f.write("@{}\n{}\n+\n{}\n".format(header, seq, qual))
        out_r.write("@{}\n{}\n+\n{}\n".format(rheader, rseq, rqual))

    ratio_dups = dups_count / items_count
    message = ("\nDereplication Complete\n\nReads/Pairs processed:\t{!s}\n"
        "Duplicates found:\t{!s} ({:.2%})\n"
        .format(items_count, dups_count, ratio_dups))
    print(message)

if __name__ == "__main__":
    main()
