#!/usr/bin/env python

"""

"""

from __future__ import print_function
import argparse
import bz2
from collections import defaultdict
import colorsys
import gzip
import pysam
import string
import sys
import zipfile

__author__ = 'Alex Hyer'
__version__ = '0.0.0a1'


class Contig:
    """Stores phylogenetic and ESOM data from short reads for a contig"""

    def __init__(self, contig_id):
        self.name = contig_id
        self.taxa_dict = defaultdict(float)
        self.chunk_numbers = []
        self.class_number = None

    def add_chunk_number(chunk_number):
        """Assign ESOM chunk numbers to contig

        :param chunk_number: NAMES file number of chunk for contig
        :type chunk_number: int or list
        """

        if type(chunk_number) is int:
            chunk_number = [chunk_number]
        self.chunk_numbers += chunk_number

    def assign_class_number(class_number):
        """Assign the contig to an ESOM class

        :param class_number: the ESOM class to assign the contig to
        :type class_number: int
        """

        self.class_number = class_number

    def add_taxa_data(taxa_name, prob_mass):
        """Add Phylosift short read data to contig information

        Note: all taxa given are assumed to be at the same taxonomic level

        :param taxa_name: The taxa the short read is associated with
        :type taxa_name: str

        :param prob_mass: The Phylosift probability mass for the read
        "type prob_mass: float
        """

        self.taxa_dict[taxa_name] += prob_mass

    def best_taxa():
        """Identify the most probable taxa for the contig and return it

        :returns: most probable taxa for the contig
        :rtype: str
        """

        try:
            taxa = max(self.taxa_dict.iteritems(), key=lambda x: x[1])[0]      
        except ValueError:
            taxa = None
        return taxa

    def possible_taxa():
        """Returns all possible taxa for contig

        :returns: all possible taxa for contig
        :rtype: view
        """

        return self.taxa_dict.keys()


def names_dict(names_handle):
    """Returns nested dictionary of NAMES file

    :returns: Dictionary as structured below
    :rtype: dict

    :param names_handle: file handle to NAMES file
    :type names_handle: File Object

    Dictionary Structure (YAML style)
    ---------------------------------

    contig_name:
        contig_chunk: chunk_number
    """

    temp_dict = defaultdict(dict)
    names_handle.readline()
    for line in names_handle:
        columns = line.strip().split('\t')
        temp_dict[columns[2]][columns[1]] = columns[0]
    return temp_dict


def rainbow_picker(scale):
    """Generates rainbow RGB values

    :returns: [scale] number of RGB tuples
    :rtype: list

    :param scale: number of RGB values to generate
    :type scale: int
    """

    hsv_tuples = [(i / scale, 1.0, 1.0) for i in range(scale)]
    rgb_tuples = map(lambda x: tuple(i * 255 for i in \
                     colorsys.hsv_to_rgb(*x)), hsv_tuples)
    return rgb_tuples


def taxa_dict(taxa_handle):
    """Returns nested dictionary of sequence_taxa_summary.txt

    :returns: Dictionary as structured below
    :rtype: dict

    :param taxa_handle: file handle to sequence_taxa_summary.txt
    :type taxa_handle: File Object

    Dictionary Structure (YAML style)
    ---------------------------------

    short_read_name:
        taxa_level: [taxa_name,probability_mass]
    """

    temp_dict = defaultdict(dict)
    taxa_handle.readline()
    for line in taxa_handle:
        columns = line.strip().split('\t')
        temp_dict[columns[0]][columns[3]] = [columns[4], columns[5]]
    return temp_dict


def x_reader(file_name):
    """Detect compression type and return appropriate file handle

    :returns: A file handle depending on file type
    :rtype: File Handle

    :param file_name: Name of file to open
    :type file_name: str

    Supports GZIP, BZIP2, and ZIP compressed files,
    returns a normal file handle if file isn't compressed.
    """

    supported_files = {
        'gz': gzip.open,
        'bz2': bz2.BZ2File,
        'zip': zipfile.ZipFile.open
    }
    last_ext = file_name.split('.')[-1]
    if last_ext in supported_files:
        return supported_files[last_ext](file_name, 'rU')
    else:
        return open(file_name, 'rU')


def main(args):
    # Instantiate each contig and assign chunk numbers
    names = names_dict(args.names)
    args.names.close()
    contigs = defaultdict(dict)
    for name in names:
        contigs[name] = Contig(name)
        chunk_numbers = [names[name][chunk] for chunk in names[name]]
        contigs[name].add_chunk_numbers(chunk_numbers)

    # Add taxonomy data to Contig based on what short reads map to them
    taxa = taxa_dict(args.taxonomy)
    args.taxonomy.close()
    unique_taxa = {'N/A': 1}
    unique_taxa_number = 2 
    for reference in args.bam.references:
        if reference in contigs:
            for read in args.bam.fetch(reference=reference):
                read_name = read[0]
                if read_name in taxa:
                    taxa_name = taxa[read_name][args.taxa_level][0]
                    prob_mass = taxa[read_name][args.taxa_level][1]
                    contigs[reference].add_taxa_data(taxa_name, prob_mass)
                    if taxa_name not in unique_taxa:
                        unique_taxa[taxa_name] = unique_taxa_number
                        unique_taxa_number += 1
    args.bam.close()

    # Assign each contig a class number based on most likely taxa
    class_file_dict = defaultdict(int)
    for contig in contigs:
        best_taxa = contigs[contig].best_taxa()
        if best_taxa is None:
            best_taxa = 'N/A'
        class_number = unique_taxa[best_taxa]
        # contigs[contig].assign_class_number(class_number)
        for chunk in contigs[contig].chunk_numbers:
            class_file_dict[chunk] = class_number

    # Color classes
    header_colors = []
    rgb_tuples = rainbow_picker(len(unique_taxa))
    for rgb_tuple in enumerate(rgb_tuples):
        color = '%{0} {1}\t{2}\t{3}'.format(rgb_tuple[0],
                                            rgb_tuple[1][0],
                                            rgb_tuple[1][1],
                                            rgb_tuple[1][2])
        header_colors.append(color)

    # Write .cls file
    taxa_output = args.output.name.replace('.cls', '.taxa')
    args.output.write('% {0}\n'.format(len(contigs)))
    args.output.write('{0}\n'.format('\n'.join(header_colors))) 
    for key in sorted(class_file_dict.keys()):
        value = class_file_dict[key]
        args.output.write('{0}\t{1}\n'.format(key, value))
    args.output.close()

    # Write .taxa file correlating class and taxonomy
    with open(taxa_output, 'w') as taxa_handle:
        taxa_handle.write('Class\tTaxonomy\n')
        for taxa in sorted(unique_taxa.items(), key=lambda x: x[1]):
            taxa_handle.write('{0}\t{1}\n'.format(taxa[1], taxa[0]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.
                                     RawDescriptionHelpFormatter)
    parser.add_argument('--bam', metavar='BAM file',
                        required=True,
                        type=pysam.AlignmentFile,
                        help='BAM file containing alignment data '
                             'for short reads that aligned to the FASTA file')
    parser.add_argument('--names', metavar='NAMES file',
                        required=True,
                        type=x_reader,
                        help='NAMES file from ESOM')
    parser.add_argument('--taxonomy',
                        required=True,
                        type=x_reader,
                        help='Phylosift sequence_taxa_summary.txt file '
                              'from Phylosift run on same short reads as '
                              'the BAM file')
    parser.add_argument('--taxa_level', metavar='Taxonomy level',
                        required=True,
                        type=string.lower,
                        help='taxonomic rank to use for color filter')
    parser.add_argument('--output', metavar='OUT file',
                        required=True,
                        type=argparse.FileType('w'),
                        help='Output file to write, ".cls" will be added')
    args = parser.parse_args()

    main(args)

    sys.exit(0)
