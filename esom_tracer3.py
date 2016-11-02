#! /usr/bin/env python

"""
Copyright:

    esom_tracer3.py Color ESOM best matches by metabolic pathways
    Copyright (C) 2016  William Brazelton, Alex Hyer

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
import bz2
from collections import defaultdict
import colorsys
import gzip
import string
import sys
import zipfile

__author__ = 'Alex Hyer'
__version__ = '0.0.1'


class Contig:
    """Stores pathway information and ESOM data for a contig"""

    def __init__(self, contig_id):
        self.chunk_numbers = []
        self.class_number = None
        self.name = contig_id
        self.pathways = []

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
        name = '_'.join(columns[2].split('_')[0:2]).strip()
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
    contigs = {}
    for name in names:
        contigs[name] = Contig(name)
        chunk_numbers = [int(names[name][chunk]) for chunk in names[name]]
        contigs[name].add_chunk_numbers(chunk_numbers)

    # Read in table per pathway
    contigs_with_pathways = 0
    contigs_not_found = 0
    unique_pathways_number = 1
    unique_pathways = {'N/A': 1}
    print('> Processing {0}'.format(args.pathways_table.name))
    for line in args.pathways_table:
        line = line.split('\t')
        pathway_contigs = line[3].split(',')[:-1]
        if line[0] not in unique_pathways:
            unique_pathways_number += 1
            unique_pathways[line[0]] = unique_pathways_number
        for contig in pathway_contigs:
            try:
                if line[0] not in contigs[contig].pathways:
                    contigs[contig].pathways.append(line[0])
                if len(contigs[contig].pathways) == 1:
                    contigs_with_pathways += 1
            except KeyError:
                contigs_not_found += 1

    print('> {0} contigs in {1} don\'t contain pathways'.format(str(contigs_not_found),
                                                                args.pathways_table.name))
    print('> {0} contigs in {1} contain pathways'.format(str(contigs_with_pathways),
                                                           args.pathways_table.name))
    print('> {0} contained {1} unique pathways'.format(args.pathways_table.name,
                                                       str(unique_pathways_number)))
    print('> Finished processing {0}'.format(args.pathways_table.name))
    args.pathways_table.close()

    # Assign each contig a class number based on pathway
    print('> Assigning classes to contigs')
    class_file_dict = defaultdict(int)
    for contig in contigs:
        if len(contigs[contig].pathways) > 1 and \
                '/'.join(contigs[contig].pathways) not in unique_pathways:
            #print('> {0} contains a unique combination of pathway assignments, '
            #      'creating new class.'.format(contig))
            unique_pathways_number += 1
            unique_pathways['/'.join(contigs[contig].pathways)] = unique_pathways_number
            contigs[contig].assign_class_number(unique_pathways[contigs[contig].pathways[0]])
            for chunk in contigs[contig].chunk_numbers:
                class_file_dict[chunk] = unique_pathways[contigs[contig].pathways[0]]
        elif contigs[contig].pathways != []:
            contigs[contig].assign_class_number(unique_pathways[contigs[contig].pathways[0]])
            for chunk in contigs[contig].chunk_numbers:
                class_file_dict[chunk] = unique_pathways[contigs[contig].pathways[0]]
        else:
            contigs[contig].assign_class_number(1)
            for chunk in contigs[contig].chunk_numbers:
                class_file_dict[chunk] = 1
    print('> {0} contigs assigned classes'.format(contigs_with_pathways))

    # Color classes
    header_colors = ['%1 255\t255\t255'] # Default class is white
    rgb_tuples = rainbow_picker(unique_pathways_number - 1) # Ignore first class
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
    pathways_output = args.output.name.replace('.cls', '.taxa')
    if not pathways_output.endswith('.taxa'):
        pathways_output += '.taxa'
    args.output.write('% {0}\n'.format(len(class_file_dict)))
    args.output.write('{0}\n'.format('\n'.join(header_colors)))
    for key in sorted(class_file_dict.keys()):
        value = class_file_dict[key]
        args.output.write('{0}\t{1}\n'.format(key, value))
    args.output.close()

    # Write .taxa file correlating class and pathways
    print('> Writing {0}'.format(pathways_output))
    with open(pathways_output, 'w') as pathway_handle:
        pathway_handle.write('Class\tTaxonomy\n')
        for pathway in sorted(unique_pathways.items(), key=lambda x: x[1]):
            pathway_handle.write('{0}\t{1}\n'.format(pathway[1], pathway[0]))
    print('> {0} unique taxa written to {1}'.format(
            str(unique_pathways_number - 1),
            pathways_output))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.
                                     RawDescriptionHelpFormatter)
    parser.add_argument('--bw',
                        action='store_true',
                        help='color all data points black and white')
    parser.add_argument('--names', metavar='NAMES file',
                        required=True,
                        type=x_reader,
                        help='NAMES file from ESOM')
    parser.add_argument('--output', metavar='OUT file',
                        required=True,
                        type=argparse.FileType('w'),
                        help='Output file to write, ".cls" will be added')
    parser.add_argument('--pathways_table', metavar='Pathways Table',
                        type=x_reader,
                        help='Pathways table from pathways2contigs')
    args = parser.parse_args()

    main(args)

    sys.exit(0)
