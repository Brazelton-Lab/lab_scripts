#! /usr/bin/env python
"""
For dereplicating sequencing reads. Can check for exact duplicates and
duplicate reads that are prefix of another. Has support for bz2 and gzip
compression.

Copyright:

    derep3.py Dereplicate exact duplicates and prefix reads
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
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
from __future__ import print_function
from __future__ import division

__author__ = "Christopher Thornton"
__date__ = "2015-12-07"
__version__ = "0.3.1"

import screed.DBConstants as DBConstants
from itertools import izip
import UserDict
import os
import io
import textwrap
import sys
import gzip
import bz2
import argparse
import hashlib
FieldTypes = (('name', DBConstants._INDEXED_TEXT_KEY),
              ('annotations', DBConstants._STANDARD_TEXT),
              ('sequence', DBConstants._STANDARD_TEXT),
              ('quality', DBConstants._STANDARD_TEXT,
               'location', DBConstants._STANDARD_TEXT))

class _screed_record_dict(UserDict.DictMixin):

    """
    Simple dict-like record interface with bag behavior.
    """

    def __init__(self, *args, **kwargs):
        self.d = dict(*args, **kwargs)

    def __getitem__(self, name):
        return self.d[name]

    def __setitem__(self, name, value):
        self.d[name] = value

    def __getattr__(self, name):
        try:
            return self.d[name]
        except KeyError:
            raise AttributeError(name)

    def __len__(self):
        return len(self.sequence)

    def keys(self):
        return self.d.keys()

def fastq_iter(handle, line=None, parse_description=True):
    """
    Iterator over the given FASTQ file handle returning records. handle
    is a handle to a file opened for reading
    """
    if line is None:
        line = handle.readline()
    line = line.strip()
    while line:
        data = _screed_record_dict()

        if not line.startswith('@'):
            raise IOError("Bad FASTQ format: no '@' at beginning of line")

        # Try to grab the name and (optional) annotations
        if parse_description:
            try:
                data['name'], data['annotations'] = line[1:].split(' ', 1)
            except ValueError:  # No optional annotations
                data['name'] = line[1:]
                data['annotations'] = ''
                pass
        else:
            data['name'] = line[1:]
            data['annotations'] = ''

        # Extract the sequence lines
        sequence = []
        data['location'] = handle.tell()
        line = handle.readline().strip()
        while not line.startswith('+') and not line.startswith('#'):
            sequence.append(line)
            line = handle.readline().strip()

        data['sequence'] = ''.join(sequence)

        # Extract the quality lines
        quality = []
        line = handle.readline().strip()
        seqlen = len(data['sequence'])
        aclen = 0
        while not line == '' and aclen < seqlen:
            quality.append(line)
            aclen += len(line)
            line = handle.readline().strip()

        data['quality'] = ''.join(quality)
        if len(data['sequence']) != len(data['quality']):
            raise IOError('sequence and quality strings must be '
                          'of equal length')

        yield data

def open_output(out_file):
    """Decide how to open an outfile for writing based on file extension"""
    extension = out_file.split('.')[-1]
    if extension == 'gz':
        return gzip.GzipFile(out_file, 'wb')
    elif extension == 'bz2':
        return bz2.BZ2File(out_file, 'wb')
    else:
        return open(out_file, 'w')

def open_file(filename):
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
    for magic, ftype in file_signatures.items():
        if file_start.startswith(magic):
            compression = ftype
            break
    if compression is 'bz2':
        sequencefile = bz2file.BZ2File(filename=bufferedfile)
    elif compression is 'gz':
        if not bufferedfile.seekable():
            raise ValueError("gziped data not streamable, pipe through zcat \
                             first")
        sequencefile = gzip.GzipFile(filename=filename)
    else:
        sequencefile = bufferedfile

    return sequencefile

def get_iterator(fin, rin, rout):
    f_iter = fastq_iter(open_file(fin))
    if rin:
        if not rout:
            print("error: argument -l/--out_reverse required with -r/--reverse")
            sys.exit(1)
        r_iter = fastq_iter(open_file(rin))
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

def get_sequence(handle, location):
    handle.seek(location)
    sequence = []
    line = handle.readline().strip()
    while not line.startswith('+') and not line.startswith('#'):
        sequence.append(line)
        line = handle.readline().strip()
    return ''.join(sequence)

def search_for_replicates(generator, forward_h, reverse_h, duplicates=None, prefix=False):
    if not duplicates:
        duplicates = []
    md5 = hashlib.md5
    records = {}
    hash_table = {}
    substring_size = 120
    for record in generator:
        paired = True if len(record) == 2 else False
        if paired:
            flocation, rlocation = (record[0].location, record[1].location)
            fseq, rseq = (record[0].sequence, record[1].sequence)
            ident = record[0].name
            # verify length of seqs is greater than substring
            key_item = fseq[:substring_size] + rseq[:substring_size]
            min_size = substring_size * 2
        else:
            flocation, rlocation = (record.location, None)
            fseq, rseq = (record.sequence, None)
            ident = record.name
            key_item = fseq[:substring_size]
            min_size = substring_size
        if len(key_item) < min_size:
            print(textwrap.fill("warning: some reads are shorter than {!s}bp. "
                "Consider applying a length filter to the dataset before "
                 "dereplication".format(substring_size), 79), file=sys.stderr)
            #sys.exit(1)
            continue
        records[ident] = (flocation, rlocation)

        hash_key = md5(key_item).hexdigest()
        if hash_key not in hash_table:
            hash_table[hash_key] = [ident]
            continue
        else:
            search_group = hash_table[hash_key]
        print(ident, hash_key, search_group)

        duplicate = None
        for search_id in search_group:
        # search item may have been deleted from if it was a duplicate
            if search_id == ident:
                sys.exit("error")
            try:
                fsearch_location, rsearch_location = records[search_id]
            except KeyError:
                continue
            fsearch = get_sequence(forward_h, fsearch_location)
            fexact, fprefix = replicate_status(fseq, fsearch, prefix)
            if fexact or fprefix:
                if rseq:
                    rsearch = get_sequence(reverse_h, rsearch_location)
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

    fname = args.in_f
    rname = args.in_r
    out_f = args.out_f
    out_r = args.out_r
    print(textwrap.fill("Starting {} with arguments: {}"
        .format(os.path.basename(__file__), ' '.join(sys.argv[1:]), 79)))

    seq_iter = get_iterator(fname, rname, out_r)
    duplicates = search_for_replicates(seq_iter, open_file(fname), open_file(rname), prefix=args.prefix)
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
