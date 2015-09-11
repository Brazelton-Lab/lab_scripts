#!/usr/bin/env python

from __future__ import print_function

"""Creates a class file highlighting 

Usage:

    esom_tracer.py [--coverage] <BAM file> <NAMES file> <OUT file>

Synopsis:

    Creates a CLASS file so that ESOM BestMatches points
    colored by whether or not they are found in the BAM file.
    By using SAM tools with the "-F 4" flag to filter the BAM file,
    then only contigs with at least one alignment will colored.

Required Arguments:

    BAM file:    A BAM file consisting of reads mapped to contigs
                 used in ESOM
    NAMES file:  The NAMES file used in ESOM for the contigs found
                 in the BAM file
    OUT file:    CLASS file to be written, importing this into ESOM
                 will color the appropriate points 
"""

import argparse
from bio_utils.file_tools.file_check import FileChecker
from collections import defaultdict
import pysam
import sys

__author__ = 'Alex Hyer'
__version__ = '1.0.0.0'


def names_dict(names_file):
    """Returns nested dictionary of NAMES file"""

    temp_dict = defaultdict(dict)
    with open(names_file, 'rU') as names_handle:
        names_handle.readline()
        for line in names_handle:
            split_line = line.strip().split('\t')
            temp_dict[split_line[2]][split_line[1]] = split_line[0]
    return temp_dict


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.
                                     RawDescriptionHelpFormatter)
    parser.add_argument('bam_file', metavar='BAM file',
                        help='BAM file containing alignment data '
                             'only for reads that aligned to a reference')
    parser.add_argument('names_file', metavar='NAMES file',
                        help='NAMES file for "upstream" ESOM map')
    parser.add_argument('fasta_file', metavar='FASTA file',
                        help='FASTA file corresponding ot NAMES file')
    parser.add_argument('out_file', metavar='OUT file',
                        help='Output file to write, ".cls" will be added')
    parser.add_argument('-c', '--coverage',
                        type=int,
                        default=50,
                        help='minimum number of non-zero base coverages')
    args = parser.parse_args()

    bamFile = FileChecker(args.bam_file)
    namesFile = FileChecker(args.names_file)
    fastaFile = FileChecker(args.fasta_file)
    outFile = FileChecker(args.out_file + '.cls')
    bamFile.read_check()
    namesFile.read_check()
    fastaFile.read_check()
    outFile.write_check()

    references = {}
    fasta_handle = pysam.FastaFile(fastaFile.name())
    with pysam.AlignmentFile(bamFile.name(), 'rb') as bam_handle:
        for reference in bam_handle.references:
            read_count = bam_handle.count(reference=reference)
            if read_count != 0:
                bases = 0
                for base in bam_handle.pileup(reference=reference):
                    bases += 1
                referenceLength = fasta_handle.get_reference_length(reference)
                baseCoverage = float(bases) / float(referenceLength)
                coverageThreshold = float(args.coverage) / 100.0
                if baseCoverage >= coverageThreshold:
                    references[reference.rsplit('_', 1)[0]] = ''
    fasta_handle.close()
    namesDict = names_dict(namesFile.name())
    classes = defaultdict(int)
    for name in namesDict:
        hit = True if name in references else False
        for subName in namesDict[name]:
            classes[int(namesDict[name][subName])] = 1 if hit is True else 0
    with open(outFile.name(), 'w') as out_handle:
        out_handle.write('% {0}\n'.format(len(references)))
        out_handle.write('%0 255\t255\t255\n')
        out_handle.write('%1 255\t0\t0\n')
        for key in sorted(classes.keys()):
            value = classes[key]
            out_handle.write('{0}\t{1}\n'.format(str(key), str(value)))

    sys.exit(0)

                                                       
