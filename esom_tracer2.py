#!/usr/bin/env python

"""

"""

from __future__ import print_function
import argparse
import bz2
from collections import defaultdict
import gzip
import pysam
from screed.openscreed import open_reader
import sys
import zipfile

__author__ = 'Alex Hyer'
__version__ = '0.0.0a1'


class Contig:
    """Stores phylogenetic data from short reads for a single contig"""

    def __init__(self, contig_id):
        self.name = contig_id
        self.taxa_dict = defaultdict(float)

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

        taxa = max(self.taxa_dict.iteritems(), key=lambda x: x[1])[0]      
        return taxa

    def possible_taxa():
        """Returns all possible taxa for contig

        :returns: all possible taxa for contig
        :rtype: view
        """

        return self.taxa_dict.keys()


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


def names_dict(names_handle):
    """Returns nested dictionary of NAMES file

    :returns: Dictionary as structured below
    :rtype: dict

    :param names_handle: file handle to NAMES file
    :type names_handle: File Object

    Dictionary Structure (YAML style)
    ---------------------------------

    contig_name:
        contig_chunk: class_number
    """

    temp_dict = defaultdict(dict)
    names_handle.readline()
    for line in names_handle:
        columns = line.strip().split('\t')
        temp_dict[columns[2]][columns[1]] = columns[0]
    return temp_dict


def main(args):
    pass


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
    parser.add_argument('--fasta', metavar='FASTA file',
                        required=True,
                        type=open_reader,
                        help='FASTA file corresponding to NAMES file')
    parser.add_argument('--taxonomy',
                        required=True,
                        type=x_reader,
                        help='Phylosift sequence_taxa_summary.txt file '
                              'from Phylosift run on same short reads as '
                              'the BAM file')
    parser.add_argument('--tax_level', metavar='Taxonomy level',
                        required=True,
                        help='taxonomic rank to use for color filter')
    parser.add_argument('--output', metavar='OUT file',
                        required=True,
                        type=argparse.FileType('w'),
                        help='Output file to write, ".cls" will be added')
    args = parser.parse_args()

    main(args)

    sys.exit(0)
