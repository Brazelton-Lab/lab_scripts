#!/usr/bin/env python

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

from __future__ import division
from __future__ import print_function
import argparse
from bio_utils.file_tools.file_check import FileChecker
from collections import defaultdict
import colorsys
import pysam
import sys

__author__ = 'Alex Hyer'
__version__ = '1.0.2.0'


def names_dict(names_file):
    """Returns nested dictionary of NAMES file"""

    temp_dict = defaultdict(dict)
    with open(names_file, 'rU') as names_handle:
        names_handle.readline()
        for line in names_handle:
            columns = line.strip().split('\t')
            temp_dict[columns[2]][columns[1]] = columns[0]
    return temp_dict


def color_taxa(names_dict, references, tax_file, tax_level):
    """Taxes names dictionary from this script and colors by taxa"""

    tax_level = tax_level.lower()
    taxa = {}
    classes = {}
    taxa_number = 1
    with open(tax_file, 'rU') as tax_handle:
        tax_handle.readline()
        for line in tax_handle:
            line = line.strip()
            columns = line.split('\t')
            if columns[0] in names_dict and 
                    columns[0] in references and
                    columns[3] == tax_level:
                if not columns[3] in taxa:
                    taxa[columns[3]] = taxa_number
                    taxa_number += 1
                classes[columns[0]] = taxa[columns[3]]
    rgb_tuples = rainbow_picker(taxa_number)
    header_colors = []
    for rgb_tuple in enumerate(rgb_tuples):
        color = '%{0} {1}\t{2}\t{3}'.format(rgb_tuple[0],
                                            rgb_tuple[1][0],
                                            rgb_tuple[1][1],
                                            rgb_tuple[1][2])
        header_colors.append(color)
    return classes, header_colors, taxa


def rainbow_picker(scale):
    """Generates rainbow RBG values"""

    hsv_tuples = [(i / scale, 1.0, 1.0) for i in range(scale)]
    rgb_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), hsv_tuples)
    return rgb_tuples


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.
                                     RawDescriptionHelpFormatter)
    parser.add_argument('--bam', metavar='BAM file',
                        help='BAM file containing alignment data '
                             'only for reads that aligned to a reference')
    parser.add_argument('--names', metavar='NAMES file',
                        help='NAMES file for "upstream" ESOM map')
    parser.add_argument('--fasta', metavar='FASTA file',
                        help='FASTA file corresponding to NAMES file')
    parser.add_argument('--taxonomy',
			help='Phylosift sequence_taxa_summary.1.txt file')
    parser.add_argument('--tax_level', metavar='Taxonomy level',
                        help='taxonomic rank to use for color filter')
    parser.add_argument('out_file', metavar='OUT file',
                        help='Output file to write, ".cls" will be added')
    parser.add_argument('-c', '--coverage',
                        type=int,
                        default=50,
                        help='minimum percentage of non-zero bases to allow')
    args = parser.parse_args()

    if args.taxonomy and not arge.tax_level:
        print('--taxonomy and --tax_level must be specified together.')
        sys.exit(1)
    references = {}
    with pysam.AlignmentFile(args.bam, 'rb') as bam_handle \
            pysam.FastaFile(args.fasta):
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
                    references[reference] = ''
    namesDict = names_dict(args.names)
    if not args.taxonomy:
        classes = defaultdict(int)
        for name in namesDict:
            hit = True if name in references else False
            for subName in namesDict[name]:
                classes[int(namesDict[name][subName])] = 1 if hit is True else 0
        header_colors = ['%0 255\t255\t255', '%1 255\t0\t0']
    else:
        classes, header_colors, taxa = color_taxa(taxonomy, references, args.taxonomy, args.tax_level)
        with open(args.out + '.taxa', 'rU') as taxa_handle:
            taxa_handle.write('Class\tTaxonomy\n')
            for key in taxa.keys():
                taxa_handle.write('{0}\t{1}\n'.format(taxa[key], key))
    with open(args.out + '.cls', 'w') as out_handle:
        out_handle.write('% {0}\n'.format(len(references)))
        out_handle.write('{0}\n'.foramt('\n'.join(header_colors)))
        for key in sorted(classes.keys()):
            value = classes[key]
            out_handle.write('{0}\t{1}\n'.format(str(key), str(value)))

    sys.exit(0)

