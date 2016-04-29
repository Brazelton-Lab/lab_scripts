#!/usr/bin/env python

"""Color ESOM data points by Phylogeny

Usage:

    esom_tracer2.py [--bam] [--names] [--taxonomy] [--taxa_level] [--output]

Synopsis:

    Takes alignment data from short reads mapped to an assembly in BAM format
    and the phylogeny of those short reads from the Phylosift
    sequence_taxa_summary.txt file to identify and color the phylogeny of each
    contig in the assembly. 
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
__version__ = '1.3.0'


class Contig:
    """Stores phylogenetic and ESOM data from short reads for a contig"""

    def __init__(self, contig_id):
        self.name = contig_id
        self.taxa_dict = defaultdict(float)
        self.chunk_numbers = []
        self.class_number = None

    def add_chunk_numbers(self, chunk_numbers):
        """Assign ESOM chunk numbers to contig

        :param chunk_number: NAMES file number of chunk for contig
        :type chunk_number: int or list of ints
        """

        if type(chunk_numbers) is int:
            chunk_numbers = [chunk_numbers]
        self.chunk_numbers += chunk_numbers

    def assign_class_number(self, class_number):
        """Assign the contig to an ESOM class

        :param class_number: the ESOM class to assign the contig to
        :type class_number: int
        """

        self.class_number = class_number

    def add_taxa_data(self, taxa_name, prob_mass):
        """Add Phylosift short read data to contig information

        Note: all taxa given are assumed to be at the same taxonomic level

        :param taxa_name: The taxa the short read is associated with
        :type taxa_name: str

        :param prob_mass: The Phylosift probability mass for the read
        "type prob_mass: float
        """

        self.taxa_dict[taxa_name] += prob_mass

    def best_taxa(self):
        """Identify the most probable taxa for the contig and return it

        :returns: most probable taxa for the contig
        :rtype: str
        """

        try:
            taxa = max(self.taxa_dict.iteritems(), key=lambda x: x[1])[0]      
        except ValueError:
            taxa = None
        return taxa

    def possible_taxa(self):
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

    Dictionary Structure (YAML format)
    ----------------------------------

    contig_name:
        contig_chunk: chunk_number
    """

    temp_dict = defaultdict(dict)
    names_handle.readline()
    for line in names_handle:
        columns = line.strip().split('\t')
        name = '-'.join(columns[2].split('_')[0:2]).strip()
        temp_dict[name][columns[1]] = columns[0]
    return temp_dict


def rainbow_picker(scale):
    """Generates rainbow RGB values

    :returns: [scale] number of RGB tuples
    :rtype: list

    :param scale: number of RGB values to generate
    :type scale: int
    """

    hsv_tuples = [(float(i) / float(scale), 1.0, 1.0) for i in range(scale)]
    rgb_tuples = map(lambda x: tuple(i * 255 for i in \
                     colorsys.hsv_to_rgb(*x)), hsv_tuples)
    return rgb_tuples


def taxa_dict(taxa_handle):
    """Returns nested dictionary of sequence_taxa_summary.txt

    :returns: Dictionary as structured below
    :rtype: dict

    :param taxa_handle: file handle to sequence_taxa_summary.txt
    :type taxa_handle: File Object

    Dictionary Structure (YAML format)
    ----------------------------------

    short_read_name:
        taxa_level: [taxa_name,probability_mass]
    """

    temp_dict = defaultdict(dict)
    taxa_handle.readline()
    for line in taxa_handle:
        columns = line.strip().split('\t')
        temp_dict[columns[0].strip().split()[0]][columns[3]] = [columns[4], columns[5]]
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
    print(' '.join(sys.argv[:]))

    # Instantiate each contig and assign chunk numbers
    print('> Processing {0}'.format(args.names.name))
    names = names_dict(args.names)
    args.names.close()
    print('> Processed {0} unique contigs from {1}'.format(str(len(names)),
                                                         args.names.name))
    contigs = defaultdict(dict)
    for name in names:
        contigs[name] = Contig(name)
        chunk_numbers = [int(names[name][chunk]) for chunk in names[name]]
        contigs[name].add_chunk_numbers(chunk_numbers)

    # Add taxonomy data to Contig based on what short reads map to them
    print('> Processing {0}'.format(args.taxonomy.name))
    taxa = taxa_dict(args.taxonomy)
    args.taxonomy.close()
    print('> Processed {0} short reads from {1}'.format(str(len(taxa)),
                                                        args.taxonomy.name))
    unique_taxa = {'N/A': 1}
    unique_taxa_number = 2 
    print('> Processing {0}'.format(args.bam.filename))
    references_match_contigs = 0
    reads_mapping_contigs = 0
    mapped_taxa_reads = 0
    for reference in args.bam.references:
        if reference in contigs:
            references_match_contigs += 1
            for read in args.bam.fetch(reference=reference):
                reads_mapping_contigs += 1
                read_name = read.query_name
                if read_name in taxa and args.taxa_level in taxa[read_name]:
                    mapped_taxa_reads += 1
                    taxa_name = taxa[read_name][args.taxa_level][0]
                    prob_mass = float(taxa[read_name][args.taxa_level][1])
                    contigs[reference].add_taxa_data(taxa_name, prob_mass)
                    if taxa_name not in unique_taxa:
                        unique_taxa[taxa_name] = unique_taxa_number
                        unique_taxa_number += 1
    args.bam.close()
    print('> {0} contigs in {1} matched contigs in {2}'.format(
            str(references_match_contigs),
            args.bam.filename,
            args.names.name))
    print('> {0} reads from {1} map to contigs in {2}'.format(
            str(reads_mapping_contigs),
            args.bam.filename,
            args.names.name))
    print('> {0} reads from {1} map to contigs in {2} and have assigned '
          'taxa from {3} at the level {4}'.format(str(mapped_taxa_reads),
                                                  args.bam.filename,
                                                  args.names.name,
                                                  args.taxonomy.name,
                                                  args.taxa_level))
    print('> Finished processing {0}'.format(args.bam.filename))

    # Assign each contig a class number based on most likely taxa
    print('> Assigning taxa to contigs')
    class_file_dict = defaultdict(int)
    contigs_with_taxa = 0
    for contig in contigs:
        best_taxa = contigs[contig].best_taxa()
        if best_taxa is None:
            best_taxa = 'N/A'
        else:
            contigs_with_taxa += 1
        class_number = unique_taxa[best_taxa]
        # contigs[contig].assign_class_number(class_number)
        for chunk in contigs[contig].chunk_numbers:
            class_file_dict[chunk] = class_number
    print('> {0} unique contigs assigned taxa'.format(str(contigs_with_taxa)))

    # Color classes
    header_colors = ['%1 255\t255\t255'] # Default class is white
    rgb_tuples = rainbow_picker(len(unique_taxa) - 1) # Ignore first class
    for rgb_tuple in enumerate(rgb_tuples):
        color = '%{0} {1}\t{2}\t{3}'.format(rgb_tuple[0] + 2, # Skip first cls
                                            int(rgb_tuple[1][0]),
                                            int(rgb_tuple[1][1]),
                                            int(rgb_tuple[1][2]))
        if args.bw:
            color = '%{0} 0\t0\t0'.format(rgb_tuple[0] + 2)
        header_colors.append(color)

    # Write .cls file
    print('> Writing {0}'.format(args.output.name))
    taxa_output = args.output.name.replace('.cls', '.taxa')
    if not taxa_output.endswith('.taxa'):
        taxa_output += '.taxa'
    args.output.write('% {0}\n'.format(len(class_file_dict)))
    args.output.write('{0}\n'.format('\n'.join(header_colors))) 
    for key in sorted(class_file_dict.keys()):
        value = class_file_dict[key]
        args.output.write('{0}\t{1}\n'.format(key, value))
    args.output.close()

    # Write .taxa file correlating class and taxonomy
    print('> Writing {0}'.format(taxa_output))
    with open(taxa_output, 'w') as taxa_handle:
        taxa_handle.write('Class\tTaxonomy\n')
        for taxa in sorted(unique_taxa.items(), key=lambda x: x[1]):
            taxa_handle.write('{0}\t{1}\n'.format(taxa[1], taxa[0]))
    print('> {0} unique taxa written to {1}'.format(
            str(unique_taxa_number - 1),
            taxa_output))

    print('> Output files written, quitting esom_tracer2.py')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.
                                     RawDescriptionHelpFormatter)
    parser.add_argument('--bam', metavar='BAM file',
                        required=True,
                        type=pysam.AlignmentFile,
                        help='BAM file containing alignment data '
                             'for short reads that aligned to the FASTA file')
    parser.add_argument('--bw',
                        action='store_true',
                        help='color all data points black and white')
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
