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
__version__ = "0.3.1"

from screed.openscreed import open_reader
from itertools import izip
import os
import textwrap
import sys
import gzip
import bz2
import argparse
import hashlib

def open_output(out_file):
    """Decide how to open an outfile for writing based on file extension"""
    extension = out_file.split('.')[-1]
    if extension == 'gz':
        return gzip.GzipFile(out_file, 'wb')
    elif extension == 'bz2':
        return bz2.BZ2File(out_file, 'wb')
    else:
        return open(out_file, 'w')
    
def get_iterator(in_f, in_r, out_r):
    f_iter = open_reader(in_f)
    if in_r:
        if not out_r:
            print("error: argument -l/--out_reverse required with -r/--reverse")
            sys.exit(1)
        r_iter = open_reader(in_r)
        return izip(f_iter, r_iter)
    else:
        return f_iter

def replicate_status(query, queried, prefix=False):
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
        if prefix:
            query_len = len(query)
            queried_len = len(queried)
            if query_len > queried_len:
                if query[:queried_len] == queried:
                    prefix_status = 2
            else:
                if query == queried[:query_len]:
                    prefix_status = 1
    return (exact_status, prefix_status)

def get_duplicate_id(query_id, search_group, seq_db, prefix=False):
    """Return header of replicate sequence"""
    fquery, rquery = seq_db[query_id]
    # check forward read first. If it is a duplicate then check reverse
    for search_id in search_group:
        # seearch item may have been deleted from if it was a duplicate
        try:
            fsearch, rsearch = seq_db[search_id]
        except KeyError:
            continue
        fexact, fprefix = replicate_status(fquery, fsearch, prefix)
        if fexact or fprefix:
            if rquery:
                rexact, rprefix = replicate_status(rquery, rsearch, prefix)
                if not rexact or rprefix:
                    # continue looking for forward duplicates (two pairs may 
                    # have prefix forward reads but different reverse reads)
                    continue
                elif ((fexact == 1) and (rexact == 1)) or \
                    ((fexact == 1) and (rprefix == 1)) or \
                    ((fprefix == 1) and (rexact == 1)) or \
                    ((fprefix == 1) and rprefix == 1):
                    return query_id
                elif ((fexact == 1) and (rprefix == 2)) or \
                    ((fprefix == 2) and (rexact == 1)) or \
                    ((fprefix == 2) and (rprefix == 2)):
                    return search_id
                # implies (fexact/fprefix == 1/2) and (rprefix/rexact== 2/1)
                else:
                    # don't know which is a prefix of which, so don't count
                    continue
            else:
                if (fexact == 1) or (fprefix == 1):
                    return query_id
                elif fprefix == 2:
                    return search_id
    return None

def search_for_replicates(generator, duplicates=None, prefix=False):
    if not duplicates:
        duplicates = []
    md5 = hashlib.md5
    records = {}
    hash_table = {}
    substring_size = 20
    for record in generator:
        paired = True if len(record) == 2 else False
        if paired:
            fseq, rseq = (record[0].sequence, record[1].sequence)
            ident = record[0].name.split()[0]
            # verify length of seqs is greater than substring
            key_item = fseq[:substring_size] + rseq[:substring_size]
            min_size = substring_size * 2
        else:
            fseq, rseq = (record.sequence, None)
            ident = record.name.split()[0]
            key_item = fseq[:substring_size]
            min_size = substring_size
        if len(key_item) < min_size:
            print(textwrap.fill("warning: some reads are shorter than {!s}bp. "
                "Consider applying a length filter to the dataset before "
                 "dereplication".format(substring_size), 79), file=sys.stderr)
            sys.exit(1)
        records[ident] = (fseq, rseq)

        hash_key = md5(key_item).digest()
        if hash_key not in hash_table:
            hash_table[hash_key] = [ident]
            continue
        else:
            search_group = hash_table[hash_key]

        duplicate = None
        for search_id in search_group:
        # search item may have been deleted from if it was a duplicate
            if search_id == ident:
                sys.exit("error")
            try:
                fsearch, rsearch = records[search_id]
            except KeyError:
                continue
            fexact, fprefix = replicate_status(fseq, fsearch, prefix)
            if fexact or fprefix:
                if rseq:
                    rexact, rprefix = replicate_status(rseq, rsearch, prefix)
                    if not rexact or rprefix:
                        # continue looking for forward dups (two pairs may have
                        # prefix forward reads but different reverse reads)
                        continue
                    elif ((fexact == 1) and (rexact == 1)) or \
                        ((fexact == 1) and (rprefix == 1)) or \
                        ((fprefix == 1) and (rexact == 1)) or \
                        ((fprefix == 1) and rprefix == 1):
                        duplicate = ident
                        break
                    elif ((fexact == 1) and (rprefix == 2)) or \
                        ((fprefix == 2) and (rexact == 1)) or \
                        ((fprefix == 2) and (rprefix == 2)):
                        duplicate = search_id
                        break
                # implies (fexact/fprefix == 1/2) and (rprefix/rexact== 2/1)
                else:
                    # don't know which is a prefix of which, so don't count
                    continue
            else:
                if (fexact == 1) or (fprefix == 1):
                    duplicate = ident
                    break
                elif fprefix == 2:
                    duplicate = search_id
                    break

        if duplicate:
            duplicates.append(duplicate)
            if duplicate != ident:
                hash_table[hash_key].remove(duplicate)
                hash_table[hash_key].append(ident)
        else:
            hash_table[hash_key].append(ident)

    return duplicates

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
    parser.add_argument('-o', '--out_forward', dest='out_f', metavar='FASTQ',
        required=True,
        type=open_output,
        help="output forward reads in fastq format [required]")
    parser.add_argument('-l', '--out_reverse', dest='out_r', metavar='FASTQ',
        type=open_output,
        help="output reverse reads in fastq format (use with -r/--reverse)")
    parser.add_argument('-p', '--prefix',
        action='store_true',
        help="replicate can be a 5' prefix of another read")
    parser.add_argument('-c', '--revcomp',
        action='store_true',
        help="replicate can be a reverse compliment of another read")
    args = parser.parse_args()

    in_f = args.in_f
    in_r = args.in_r
    out_f = args.out_f
    out_r = args.out_r
    print(textwrap.fill("Starting {} with arguments: {}"
        .format(os.path.basename(__file__), ' '.join(sys.argv[1:]), 79)))

    seq_iter = get_iterator(in_f, in_r, out_r)
    duplicates = search_for_replicates(seq_iter, prefix=args.prefix)
    dups_count = len(duplicates)
    print("Num duplicates: {!s}".format(dups_count))

    # time to write the output
    sys.exit(0)
    seq_iter = get_iterator(in_f, in_r, out_r)
    items_count = 0
    for record in seq_iter:
        items_count += 1
        if in_r:
            seq, rseq = (record[0].sequence, record[1].sequence)
            qual, rqual = (record[0].quality, record[1].quality)
            header, rheader = (record[0].name, record[1].name)
        else:
            seq = record.sequence
            qual = record.quality
            header = record.name

        ident = header.split()[0]
        if ident in duplicates:
            continue
        out_f.write("@{}\n{}\n+\n{}\n".format(header, seq, qual))
        if out_r:
            out_r.write("@{}\n{}\n+\n{}\n".format(rheader, rseq, rqual))

    dups_count = len(duplicates)
    ratio_dups = dups_count / items_count
    print("\nDereplication Complete\n\nReads/Pairs processed:\t{!s}\nDuplicates "
        "found:\t{!s} ({:.2%})\n\n".format(items_count, dups_count, ratio_dups))

if __name__ == "__main__":
    main()
