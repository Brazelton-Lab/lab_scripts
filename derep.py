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
import argparse
import hashlib

def open_output(out_file):
    """decide how to open an outfile for writing based on file extension"""
    extension = out_file.split('.')[-1]
    if extension == 'gz':
        return gzip.open(out_file, 'wb')
    else:
        return open(out_file, 'w')

def replicate_status(query, search, prefix=False):
    status = 0
    if prefix:
        if query[:len(search)] == search:
            status = 2
        elif search[:len(query)] == query:
            status = 1
    else:
        if query == search:
            status = 1
    return status

def search_replicates(record, key, unique_db, prefix=False):
    # check forward read first. If it is a duplicate then check reverse
    header = None

    fquery = record[0]
    query_ident = record[2]
    for i, (fsearch, rsearch) in enumerate(unique_db[key]['seq']):
        fstatus = replicate_status(fquery, fsearch, prefix)
        if fstatus:
            rquery = record[1]
            if rquery:
                rstatus = replicate_status(rquery, rsearch, prefix)
                if not rstatus:
                    # continue looking (two pairs my have prefix forward reads but different reverse reads)
                    continue
                elif (fstatus == 1) and (rstatus == 1):
                    return query_ident
                elif (fstatus == 2) and (rstatus == 2):
                    search_ident = unique_db[key]['ids'][i]
                    return search_ident
                # implied ((fstatus == 1) and (rstatus == 2)) or ((fstatus == 2) and (rstatus == 1))
                else:
                    # return largest of combined
                    search_ident = unique_db[key]['ids'][i]
                    query_length = len(fquery + rquery)
                    search_length = len(fsearch + rsearch)
                    if query_length > search_length:
                        return query_ident
                    else:
                        search_ident = unique_db[key]['ids'][i]
                        return search_ident
            else:
                if fstatus == 1:
                    return query_ident
                elif fstatus == 2:
                    search_ident = unique_db[key]['ids'][i]
                    return search_ident
    return header

def main():
    parser = argparse.ArgumentParser(description="Remove exact or prefix duplicates from fastq files")
    parser.add_argument('-f', '--forward', dest='in_f',
        required=True,
        type=open_reader,
        help="paired forward reads in fastq format. Will treat as single-end "
            "reads if the pair is not provided with argument -r/--reverse")
    parser.add_argument('-r', '--reverse', dest='in_r',
        type=open_reader,
        help="paired reverse reads in fastq format")
    parser.add_argument('-o', '--out_forward', dest='out_f',
        required=True,
        type=open_output,
        help="output forward reads in fastq format")
    parser.add_argument('-l', '--out_reverse', dest='out_r',
        type=open_output,
        help="output reverse reads in fastq format (if applicable)")
    parser.add_argument('-p', '--prefix',
        action='store_true',
        help="replicate can be a 5' prefix of another read")
    args = parser.parse_args()

    md5 = hashlib.md5
    if args.in_r:
        if not args.out_r:
            print("error: argument -l/--out_reverse required with -r/--reverse")
            sys.exit(1)
        else:
            seq_iter = izip(args.in_f, args.in_r)
    else:
        seq_iter = args.in_f

    items_count = 0
    duplicates = []
    uniques = {}
    for record in seq_iter:
        substring = 10
        if args.in_r:
            f_seq, r_seq = (record[0].sequence, record[1].sequence)
            # verify length of seqs is greater than substring
            ident = record[0].name
            hash_key = md5(f_seq[:substring] + r_seq[:substring]).hexdigest()
        else:
            f_seq, r_seq = (record.sequence, None)
            ident = record.name
            hash_key = md5(f_seq[:substring]).hexdigest()

        if hash_key not in uniques:
            uniques[hash_key] = {'seq': [(f_seq, r_seq)], 'ids': [ident]}
            items_count += 1
            continue

        query = (f_seq, r_seq, ident)
        duplicate = search_replicates(query, hash_key, uniques, args.prefix)
        if duplicate:
            duplicates.append(duplicate)    
            if duplicate != ident:
                # find entry for search and replace with query
                try:
                    position = uniques[hash_key]['ids'].index(duplicate)
                except:
                    print("unknown error in database")
                    sys.exit(1)
                uniques[hash_key]['ids'][position] = duplicate
                uniques[hash_key]['seq'][position] = (f_seq, r_seq)
        else:
            uniques[hash_key]['seq'].append((f_seq, r_seq))
            uniques[hash_key]['ids'].append(ident)

        items_count += 1

    # unfortunately have to read through files again
    if args.in_r:
        for f_record, r_record in izip(args.in_f, args.in_r):
            if f_record.name not in duplicates:
                args.out_f.write("@{} {}\n{}\n+\n{}\n".format(f_record.name, 
                    f_record.annotations, f_record.sequence, f_record.quality))
                args.out_r.write("@{} {}\n{}\n+\n{}\n".format(r_record.name, 
                    r_record.annotations, r_record.sequence, r_record.quality))
        args.out_r.close()
        args.in_r.close()
    else:
        for record in args.in_f:
            if record.name not in duplicates:
                args.out_f.write("@{} {}\n{}\n+\n{}\n".format(f_record.name, 
                    f_record.annotations, f_record.sequence, f_record.quality))
    args.out_f.close()
    args.in_f.close()

    dups_count = len(duplicates)
    ratio_dups = dups_count / items_count
    print("\nDereplication Complete\n\nReads/Pairs processed:\t{!s}\nDuplicates "
        "found:\t{!s} ({:.2%})\n".format(items_count, dups_count, ratio_dups))
    print("unique hash: {}".format(str(unique_hash)))

if __name__ == "__main__":
    main()
