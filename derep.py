#! /usr/bin/env python
"""
For dereplicating sequencing reads. Can check for exact duplicates and 
duplicate reads that are prefix of another. Has support for bz2 and gzip 
compression.
"""
from __future__ import print_function
from __future__ import division

from screed.openscreed import open_reader
from itertools import izip
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

def search_replicates(query_id, key, seq_db, prefix=False):
    """Return header of replicate if one is found"""
    fquery, rquery = seq_db[key][query_id]['sequence']
    # check forward read first. If it is a duplicate then check reverse
    for search_id in seq_db[key]:
        if query_id == search_id:
            continue
        # don't bother checking known duplicates
        if not seq_db[key][search_id]['unique']:
            continue
        fsearch, rsearch = seq_db[key][search_id]['sequence']
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
                    print("query: {}, search: {}".format(query_id, search_id))
                    continue
            else:
                if (fexact == 1) or (fprefix == 1):
                    return query_id
                elif fprefix == 2:
                    return search_id
    return None

def main():
    parser = argparse.ArgumentParser(description="Remove exact or prefix duplicates from fastq files")
    parser.add_argument('-f', '--forward', dest='in_f',
        required=True,
        type=str,
        help="paired forward reads in fastq format. Will treat as single-end "
            "reads if the pair is not provided with argument -r/--reverse")
    parser.add_argument('-r', '--reverse', dest='in_r',
        type=str,
        help="paired reverse reads in fastq format")
    parser.add_argument('-o', '--out_forward', dest='out_f',
        required=True,
        type=open_output,
        help="output forward reads in fastq format")
    parser.add_argument('-l', '--out_reverse', dest='out_r',
        type=open_output,
        help="output reverse reads in fastq format (use with -r/--reverse)")
    parser.add_argument('-p', '--prefix',
        action='store_true',
        help="replicate can be a 5' prefix of another read")
    parser.add_argument('-c', '--revcomp',
        action='store_true',
        help="replicate can be a reverse compliment of another read")
    args = parser.parse_args()

    md5 = hashlib.md5
    in_f = args.in_f
    in_r = args.in_r
    out_f = args.out_f
    out_r = args.out_r
    substring_size = 20

    seq_iter = get_iterator(in_f, in_r, out_r)

    records = {}
    items_count = 0
    for record in seq_iter:
        if args.in_r:
            fseq, rseq = (record[0].sequence, record[1].sequence)
            ident = record[0].name.split()[0]
            # verify length of seqs is greater than substring
            key_item = fseq[:substring_size] + rseq[:substring_size]
        else:
            fseq, rseq = (record.sequence, None)
            ident = record.name.split()[0]
            key_item = fseq[:substring_size]

        hash_key = md5(key_item).digest()
        if hash_key not in records:
            records[hash_key] = {ident: {'sequence': (fseq, rseq),
                'unique': True}}
            items_count += 1
            continue
        else:
            records[hash_key][ident] = {'sequence': (fseq, rseq),
                'unique': True}

        duplicate = search_replicates(ident, hash_key, records, args.prefix)
        if duplicate:
            records[hash_key][duplicate]['unique'] = False

        items_count += 1

    seq_iter = get_iterator(in_f, in_r, out_r)
    uniques_count = 0
    for record in seq_iter:
        if in_r:
            seq, rseq = (record[0].sequence, record[1].sequence)
            qual, rqual = (record[0].quality, record[1].quality)
            header, rheader = (record[0].name, record[1].name)
            key_item = seq[:substring_size] + rseq[:substring_size]
        else:
            seq = record.sequence
            qual = record.quality
            header = record.name
            key_item = seq[:substring_size]

        ident = header.split()[0]
        key = md5(key_item).digest()
        # only print the uniques
        if not records[key][ident]['unique']:
                continue
        out_f.write("@{}\n{}\n+\n{}\n".format(header, seq, qual))
        if out_r:
            out_r.write("@{}\n{}\n+\n{}\n".format(rheader, rseq, rqual))

        uniques_count += 1

    dups_count = items_count - uniques_count
    ratio_dups = dups_count / items_count
    print("\nDereplication Complete\n\nReads/Pairs processed:\t{!s}\nDuplicates "
        "found:\t{!s} ({:.2%})\n".format(items_count, dups_count, ratio_dups))

if __name__ == "__main__":
    main()

