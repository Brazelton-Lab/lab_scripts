#! /usr/bin/env python

"""
Screen the results of a homology search using bit-score thresholds, 
alternative phenotype-conferring snps, or other scoring metrics.

Usage:
    homology_screen [options] [-o out.csv] in.csv

    homology_screen [options] [-o out.csv.gz] in.csv.gz

Required input is a tabular file of pairwise alignments (B6 format). Optional 
inputs directing how pairwise alignments should be screened are a tabular SNP 
file containing alternative phenotypes along with the location and nature of 
the SNP and/or a relational database in JSON format containing an appropriate 
scoring threshold (e.g. bitscore thresholds). The compression algorithm is 
automatically detected for input files based on file extension. To compress 
the output, add the appropriate file extensions to the output file name 
(e.g. .gz, .bz2). Use /dev/stdin to indicate that input is from standard 
input. Similarly, leaving out '--out' will result in the output being 
sent to standard output (stdout).

Copyright:

    homolog_screen  Screen the results of a homology search using bit-score thresholds, alternative phenotype-conferring snps, or other scoring metrics.

    Copyright (C) 2016  William Brazelton <comma-separated list of authors>

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

import argparse
from bz2 import BZ2File
from gzip import GzipFile
import json
import os
import sys
import textwrap
from time import time

__author__ = "Christopher Thornton"
__license__ = 'GPLv2'
__maintainer__ = 'Christopher Thornton'
__status__ = "Alpha"
__version__ = "0.0.1"


class FormatError(Exception):
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)


class Open(argparse.Action):
    """Argparse Action that detects and opens compressed files for rw

    Attributes:
        option_strings (list): list of str giving command line flags that
                               call this action

        dest (str): Namespace reference to value

        mode (str): mode to pass to (de)compression algorithm

        nargs (bool): True if multiple arguments specified

        **kwargs (various): optional arguments to pass to argparse and algo
    """

    def __init__(self, option_strings, dest, mode='rb', nargs=None, **kwargs):
        """Initialize class and spawn self as Base Class w/o nargs

        Warns:
            ImportError: if Open cannot import a compression library,
                         it warns the user that it cannot open the
                         corresponding file type

        Raises:
            ValueError: if nargs is not None, Open does not accept nargs
        """
        from importlib import import_module
        from warnings import warn
        from os import linesep
        # Only accept a single value to analyze
        if nargs is not None:
           raise ValueError('nargs not allowed for Open')

        # Call self again but without nargs
        super(Open, self).__init__(option_strings, dest, **kwargs)

        # Store and establish variables used in __call__
        self.kwargs = kwargs
        self.mode = mode.lower().strip()
        self.modules = {}

        modules_to_import = {
            'bz2': 'BZ2File',
            'gzip': 'GzipFile',
            'lzma': 'LZMAFile'
        }

        # Dynamically import compression libraries and warn about failures
        for mod, _class in modules_to_import.items():
            try:
                self.modules[_class] = getattr(import_module(mod), _class)
            except (ImportError, AttributeError) as e:
                self.modules[_class] = open
                warn('Cannot process {0} files due to following error:'
                     '{1}{2}{1}You will need to install the {0} library to '
                     'properly use these files. Currently, such files will '
                     'open in text mode.'.format(mod, linesep, e))
    # Credits: https://stackoverflow.com/questions/13044562/
    # python-mechanism-to-identify-compressed-file-type-and-uncompress
    def __call__(self, parser, namespace, value, option_string=None, **kwargs):
        """Detects and opens compressed files
        Args:
            parser (ArgumentParser): parser used to generate values

            namespace (Namespace): namespace to set values for

            value (str): actual value specified by user

            option_string (str): argument flag used to call this function

            **kwargs (various): optional arguments later passed to the
                                compression algorithm
        """
        from inspect import getfullargspec
        import io

        filename = value  # For readability

        algo = io.open  # Default to plaintext

        algo_map = {
            'bz2': self.modules['BZ2File'],
            'gz':  self.modules['GzipFile'],
            'xz':  self.modules['LZMAFile']
        }

        # Base compression algorithm on file extension
        ext = value.split('.')[-1]
        try:
            algo = algo_map[ext]
        except KeyError:
            pass

        # Filter all **kwargs by the args accepted by the compression algo
        algo_args = set(getfullargspec(algo).args)
        good_args = set(self.kwargs.keys()).intersection(algo_args)
        _kwargs = {arg: self.kwargs[arg] for arg in good_args}


        # Open the file using parameters defined above and store in namespace
        try:
            handle = algo(value, mode=self.mode, **_kwargs)
        except ValueError:
            mode = self.mode.lstrip('U')[0]
            handle = io.TextIOWrapper(algo(value, mode=mode, **_kwargs), encoding='utf-8')

        setattr(namespace, self.dest, handle)


class B6Entry:
    """A simple class to store data from B6/M8 entries and write them


    Attributes:
        query (str): query ID (sequence aligned with)

        subject (str): subject ID (sequence aligned to)

        perc_identical (float): percent of query and subject sequences that are
            identical

        align_len (int): length of alignment

        mismatches (int): number of mismatches in alignment

        gaps (int): number of gaps in alignment

        query_start (int): alignment start position in query sequence

        query_end (int): alignment end position in query sequence

        subject_start (int): alignment start position in subject sequence

        subject_end (int): alignment end position in subject sequence

        evalue (float): E-value of alignment

        bit_score (float): Bit score of alignment

        others (dict): a dictionary to store values of custom format specifiers
    """
    def __init__(self):
        """Initialize variables to store B6/M8 entry data"""

        self.query = None
        self.subject = None
        self.perc_identical = None
        self.align_len = None
        self.mismatches = None
        self.gaps = None
        self.query_start = None
        self.query_end = None
        self.subject_start = None
        self.subject_end = None
        self.evalue = None
        self._evalue_str = None  #store original formatting of E-value
        self.bit_score = None
        self.others = None  #store dictionary of additional format specifiers

    def write(self):
        """Return B6/M8 formatted string

        Returns:
            str: B6/M8 formatted string containing entire B6/M8 entry
        """

        if self.others:
            other_specs = "\t{!s}".format('\t'.join([self.others[i] for i in \
                          sorted(self.others)]))
        else:
            other_specs = ''

        return '{0}\t{1}\t{2}\t{3}\t{4}\t' \
               '{5}\t{6}\t{7}\t{8}\t{9}\t' \
               '{10}\t{11}{12}{13}'.format(self.query,
                                           self.subject,
                                           str(self.perc_identical),
                                           str(self.align_len),
                                           str(self.mismatches),
                                           str(self.gaps),
                                           str(self.query_start),
                                           str(self.query_end),
                                           str(self.subject_start),
                                           str(self.subject_end),
                                           self._evalue_str,
                                           other_specs,
                                           str(self.bit_score),
                                           os.linesep)


def b6_iter(handle, start_line=None, out_header=None):
    """Iterate over B6/M8 file and return B6/M8 entries

    Args:
        handle (file): B6/M8 file handle, can be any iterator so long as it
            it returns subsequent "lines" of a B6/M8 entry

        start_line (str): Next B6/M8 entry, if 'handle' has been partially read
            and you want to start iterating at the next entry, read the next
            B6/M8 entry and pass it to this variable when  calling b6_iter.
            See 'Examples.'
    Yields:
        B6Entry: class containing all B6/M8 data

    Examples:
        The following two examples demonstrate how to use b6_iter.
        Note: These doctests will not pass, examples are only in doctest
        format as per convention. bio_utils uses pytests for testing.

        >>> for entry in b6_iter(open('test.b6out')):
        ...     print(entry.query)  # Print Query ID
        ...     print(entry.subject)  # Print Subject ID
        ...     print(entry.perc_identical)  # Print % identity between seqs
        ...     print(entry.mismatches)  # Print number of mismathces in align
        ...     print(entry.gaps)  # Print number of gaps in alignment
        ...     print(entry.query_start)  # Print start of alignment on query
        ...     print(entry.query_end)  # Print end of alignment on query
        ...     print(entry.subject_start)  # Print start of align on subject
        ...     print(entry.subject_end)  # Print end of alignment on subject
        ...     print(entry.evalue)  # Print E-value of alignment
        ...     print(entry.bit_score)  # Print Bit score of alignment
        ...     print(entry.write())  # Print entry B6 entry

        >>> b6_handle = open('test.b6out')
        >>> next(b6_handle)  # Skip first line/entry
        >>> next_line = next(b6_handle)  # Store next entry
        >>> for entry in b6_iter(b6_handle, start_line=next_line):
        ...     print(entry.query)  # Print Query ID
        ...     print(entry.subject)  # Print Subject ID
        ...     print(entry.perc_identical)  # Print % identity between seqs
        ...     print(entry.mismatches)  # Print number of mismathces in align
        ...     print(entry.gaps)  # Print number of gaps in alignment
        ...     print(entry.query_start)  # Print start of alignment on query
        ...     print(entry.query_end)  # Print end of alignment on query
        ...     print(entry.subject_start)  # Print start of align on subject
        ...     print(entry.subject_end)  # Print end of alignment on subject
        ...     print(entry.evalue)  # Print E-value of alignment
        ...     print(entry.bit_score)  # Print Bit score of alignment
        ...     print(entry.write())  # Print entry B6 entry
    """

    # Speed tricks: reduces function calls
    split = str.split
    strip = str.strip

    next_line = next

    if start_line is None:
        line = next(handle)  # Read first B6/M8 entry
    else:
        line = start_line  # Set header to given header

    # Check if input is text or bytestream
    if (isinstance(line, bytes)):
        def next_line(i):
            return next(i).decode('utf-8')

        line = strip(line.decode('utf-8'))
    else:
        line = strip(line)

    # Required format specifiers
    required_specs = ['qend', 'mismatch', 'pident', 'qaccver', 'qstart', \
                      'sstart', 'bitscore', 'evalue', 'gapopen', 'send', \
                      'length', 'saccver']

    # Check if first line is header
    if line.startswith('#'):  #file in custom format
        split_line = split(line[1:], '\t')

        header = {}
        for spec in split_line:
            header[spec] = split_line.index(spec)

        # verify that required keys are all there
        for spec in required_specs:
            if spec not in header:
                print(header, file=sys.stderr)
                raise FormatError("Required format specifier '{}' is "\
                                  "missing from the header".format(spec))

        line = strip(next_line(handle))
    else:  #file in default format
        header = {
                  'qaccver': 0, 
                  'saccver': 1, 
                  'pident': 2, 
                  'length': 3, 
                  'mismatch': 4, 
                  'gapopen': 5, 
                  'qstart': 6, 
                  'qend': 7, 
                  'sstart': 8, 
                  'send': 9,
                  'evalue': 10,
                  'bitscore': 11
                 }

    if out_header:
        other_specs = '\t'.join(sorted([i for i in header if i not in \
                                required_specs]))

        if other_specs:
            out_header.write("#{}\t{}".format('\t'.join(required_specs), \
                             other_specs))
        else:
            out_header.write("#{}".format('\t'.join(required_specs)))


    # A manual 'for' loop isn't needed to read the file properly and quickly,
    # unlike fasta_iter and fastq_iter, but it is necessary begin iterating
    # partway through a file when the user gives a starting line.
    try:  # Manually construct a for loop to improve speed by using 'next'

        while True:  # Loop until StopIteration Exception raised

            split_line = split(line, '\t')

            data = B6Entry()
            data.query = split_line[header['qaccver']]
            data.subject = split_line[header['saccver']]
            data.perc_identical = float(split_line[header['pident']])
            data.align_len = int(split_line[header['length']])
            data.mismatches = int(split_line[header['mismatch']])
            data.gaps = int(split_line[header['gapopen']])
            data.query_start = int(split_line[header['qstart']])
            data.query_end = int(split_line[header['qend']])
            data.subject_start = int(split_line[header['sstart']])
            data.subject_end = int(split_line[header['send']])
            data.evalue = float(split_line[header['evalue']])
            data._evalue_str = split_line[header['evalue']]
            data.bit_score = float(split_line[header['bitscore']])

            # Add additional format specifiers as a dictionary if custom format
            other_specs = [i for i in header if i not in required_specs]

            others = {}
            for other_spec in other_specs:
                others[other_spec] = split_line[header[other_spec]]

            data.others = others

            line = strip(next_line(handle))  # Raises StopIteration at EOF

            yield data

    except StopIteration:  # Yield last B6/M8 entry
        yield data


def parse_snps_file(in_h, line=None):
    """
    Parse an SNPs file. Each line should have four columns, corresponding to
    accession number, position of SNP, wildtype base composition, and mutant 
    base composition. An opional header is allowed.
    """

    append = list.append
    join = str.join
    strip = str.strip

    if line is None:
        line = next(in_h)  # Read header

    # Check if input is text or bytestream
    if (isinstance(line, bytes)):
        def next_line(i):
            return next(i).decode('utf-8')

        line = line.decode('utf-8')
    else:
        next_line = next

    # Check if first line is header
    if line.startswith('#'):
        split_line = strip(next_line(in_h)).split('\t')
    else:
        split_line = strip(line).split('\t')

    snp_dict = {}
    try:

        while True:  #loop until StopIteration raised

            try:
                acc, mt, pt, position, wildtype, mutant = split_line
            except ValueError:
                raise FormatError("{} requires a tabular SNP file with six "
                                  "columns, corresponding to accession #, "
                                  "model type, parameter type, position, and "
                                  "wildtype and mutatant alleles"\
                                  .format(os.path.basename(__file__)))

            try:
                snp_dict[acc].append((position, wildtype, mutant))
            except KeyError:
                snp_dict[acc] = [(position, wildtype, mutant)]

            split_line = strip(next_line(in_h)).split('\t')

    except StopIteration:

        return snp_dict


def screen_alignment_quality(hit, e=10, ident=0, cov=0, score=0):
    """
    Screen hit based on quality metrics of the alignment.
    """
    condition = (hit.evalue <= e and hit.perc_identical >= ident and \
                 hit.align_len >= cov and hit.bit_score >= float(score))

    if condition:
        return True
    else:
        return False


def screen_snp(hit, aro, snps):
    """
    Screen hit for alternative phenotype-conferring SNPs using a dictionary of
    accession number, SNP value pairs.
    """
    try:
        mutations = snps[aro]
    except KeyError:
        print("{} not found in snp database".format(aro), file=sys.stderr)

        return True
    else:
        try:
            qseq = hit.others['qseq']
        except KeyError:
            print("error: the format specifier 'qseq' is required for "
                  "secondary screening of SNPs".format(), file=sys.stderr)
            sys.exit(1)

        try:
            sseq = hit.others['sseq']
        except KeyError:
            print("error: the format specifier 'sseq' is required for "
                  "secondary screening of SNPs".format(), file=sys.stderr)
            sys.exit(1)

        start = int(hit.subject_start)
        end = int(hit.subject_end)

        for mutation in mutations:
            pos, wt, sub = mutation
            pos = int(pos)

            if pos not in list(range(start, end + 1)):
                return False

            else:
                rel_pos = pos - start  #position of snp relative to alignment

                # Handle any gaps
                i = 0
                while i != (rel_pos + 1):
                    if sseq[i] == '-':
                        rel_pos += 1
                        i += 1
                    else:
                        i += 1

                # Handle database inconsistencies
                if sseq[rel_pos] != wt:
                    print("warning: possible inconsistency in database: "
                          "subject {} with residue {} at SNP position {!s} "
                          "does not match wildtype residue {}"\
                          .format(hit.subject, sseq[rel_pos], pos, wt), \
                          file=sys.stderr)
                    return False

                if qseq[rel_pos] == sub:
                    return True
                else:
                    return False

def program_info(prog, args, version):
    print("{} {!s}".format(prog, version), file=sys.stderr)
    print(textwrap.fill("Command line parameters: {}".format(' '.join(args)), 79), file=sys.stderr)
    print("", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('inhandle',
        metavar='in.b6',
        type=str,
        action=Open,
        mode='rb',
        default=sys.stdin,
        help="input best-hit alignments in B6/M8 format. A header, "
             "differentiated by a starting # character, is required whenever "
             "format specifiers other than the default BLAST+ format "
             "specifiers are used")
    parser.add_argument('-o', '--out',
        metavar='out.b6',
        type=str,
        action=Open,
        mode='wt',
        default=sys.stdout,
        help="output best-hit alignments in B6/M8 format passing screening "
             "[default: output to stdout]")
    screen_method = parser.add_argument_group(title="screening parameters")
    screen_method.add_argument('-e', '--evalue',
        type=float,
        default=10,
        help="maximum expect value allowed for a hit to be retained [default: "
             "10]")
    screen_method.add_argument('-m', '--meta',
        metavar='meta.json',
        type=str,
        action=Open,
        mode='rb',
        help="gene metadata file containing bit-score or other alignment "
             "scoring thresholds")
    screen_method.add_argument('-v', '--value',
        type=str,
        help="entry value in the gene metadata file to use as the alignment "
             "quality threshold")
    screen_method.add_argument('-i', '--identity',
        type=float,
        default=0,
        help="minimum percent identity required to retain a hit [default: 0]")
    screen_method.add_argument('-l', '--length',
        dest='aln_len',
        type=int,
        default=0,
        help="minimum alignment length required to retain a hit [default: 0]. "
             "Should ideally be used with --identity, as alignment length is "
             "not useful on its own")
    screen_method.add_argument('-s', '--snps',
        metavar='snps.csv',
        type=str,
        action=Open,
        mode='rb',
        help="tabular file containing SNPs (position, along with the wildtype "
             "and mutant alleles) that confer alternative phenotypes to an "
             "organism. This step can only be used with -m/--meta and will be "
             "performed last if any additional screening for alignment "
             "quality is also specified")
    parser.add_argument('--version',
        action='version',
        version='%(prog)s ' + __version__)
    args = parser.parse_args()
    all_args = sys.argv[1:]

    program_info(os.path.basename(__file__), all_args, __version__)

    if not (args.meta or args.evalue or args.snps or args.aln_len or \
            args.identity):
        parser.error("error: one or more of the following arguments are "
                     "required: -s/--snps, -e/--evalue, -m/--meta, "
                     "-i/--identity, or -l/--length")

    if args.snps and not args.meta:
        parser.error("error: -m/--meta required when -s/--snps supplied")

    if args.value and not args.meta:
        parser.error("error: -m/--meta required when -v/--value supplied")

    # Track program run-time
    start_time = time()


    # Assign variables based on arguments supplied by the user
    if args.meta:
        meta_data = json.load(args.meta)
    else:
        meta_data=None

    snps = parse_snps_file(args.snps) if args.snps else None


    # Screen hits for alignment quality and/or mutant alleles
    for processed_total, hit in enumerate(b6_iter(args.inhandle)):

        if meta_data:
            try:
                sub_entry = meta_data[hit.subject]
            except KeyError:
                score = 0
            else:
                if args.value:
                    try:
                        score = sub_entry[args.value]
                    except KeyError:
                        print("error: entry value {} not found in gene meta data "
                              "file".format(args.value), file=sys.stderr)
                        sys.exit(1)
                else:
                    score = 0
        else:
            score = 0

        q_pass = screen_alignment_quality(hit, e=args.evalue, score=score, ident=args.identity, cov=args.aln_len)

        if q_pass:
            if snps:
                try:
                    aro_acc = sub_entry['gene_family']
                except UnboundLocalError:
                    print("No entry in gene metadata file for {}"\
                          .format(hit.subject), file=sys.stderr)
                    sys.exit(1)

                s_pass = screen_snp(hit, aro_acc, snps)

                if s_pass:
                    args.out.write(hit.write())

            else:
                args.out.write(hit.write())


    # Calculate and print program run-time
    end_time = time()
    total_time = (end_time - start_time) / 60.0
    print("It took {:.2e} minutes to process {!s} records\n"\
          .format(total_time, processed_total), file=sys.stderr)


if __name__ == "__main__":
    main()
    sys.exit(0)

