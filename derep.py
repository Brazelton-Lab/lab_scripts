#! /usr/bin/env python
"""dereplication program using paired-end reads"""

from __future__ import print_function
from __future__ import division
from screed.openscreed import open_reader
from itertools import izip
import gzip
import argparse
import hashlib

def open_output(out_file):
    extension = out_file.split('.')[-1]
    if extension == 'gz':
        return gzip.open(out_file, 'wb')
    else:
        return open(out_file, 'wb')


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-f', '--forward',
        required=True,
        type=open_reader,
        help="paired forward reads in fastq format")
    parser.add_argument('-r', '--reverse',
        required=True,
        type=open_reader,
        help="paired reverse reads in fastq format")
    parser.add_argument('-o', '--out_forward',
        required=True,
        type=open_output,
        help="output forward reads in fastq format")
    parser.add_argument('-l', '--out_reverse',
        required=True,
        type=open_output,
        help="output reverse reads in fastq format")
    args = parser.parse_args()

    md5 = hashlib.md5
    in_f = args.forward
    in_r = args.reverse
    out_f = args.out_forward
    out_r = args.out_reverse

    unique_seqs = {}
    num_pairs = 0
    num_dups = 0

    for f_record, r_record in izip(in_f, in_r):
        f_seq, r_seq = (f_record.sequence, r_record.sequence)
        tetra = (f_seq[0:3], r_seq[0:3])
        unique_item = md5(f_seq + r_seq).digest()
        unique = False

        if tetra not in unique_seqs:
            unique_seqs[tetra] = [unique_item]
            unique = True
        else:
            if unique_item not in unique_seqs[tetra]:
                unique_seqs[tetra].append(unique_item)
                unique = True

        if unique:
            out_f.write("@{} {}\n{}\n+\n{}\n".format(f_record.name, 
                f_record.annotations, f_record.sequence, f_record.quality))
            out_r.write("@{} {}\n{}\n+\n{}\n".format(r_record.name, 
                r_record.annotations, r_record.sequence, r_record.quality))
        else:
            num_dups += 1

        num_pairs += 1

    ratio_dups = num_dups / num_pairs
    print("\nDereplication Complete\n\nRead pairs processed:\t{!s}\nDuplicates "
        "found:\t{!s} ({:.2%})\n".format(num_pairs, num_dups, ratio_dups))

    out_f.close()
    out_r.close()

if __name__ == "__main__":
    main()
