#! /usr/bin/env python

"""Make a class file using bins from anvi_converter

Copyright:

    esom_tracer4.py Color ESOM best matches by bin

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

import argparse
import bz2
from collections import defaultdict
import colorsys
import gzip
import sys
import zipfile

__author__ = 'Alex Hyer'
__email__ = 'theonehyer@gmail.com'
__license__ = 'GPLv3'
__maintainer__ = 'Alex Hyer'
__status__ = 'Production'
__version__ = '1.0.1'


class Contig:
    """Stores pathway information and ESOM data for a contig"""

    def __init__(self, contig_id):
        self.chunk_numbers = []
        self.class_number = None
        self.name = contig_id

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
        name = '_'.join(columns[1].split('_')[0:2]).strip()
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
    rgb_tuples = map(lambda x: tuple(i * 255 for i in
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
        contigs[name].assign_class_number(0)  # Default class

    classes = [0]
    count = 0
    print('> Processing {0}'.format(args.bins.name))
    for line in args.bins:
        line = line.strip()
        line = line.split('\t')
        if line[1] not in classes:
            classes.append(line[1])
        if line[0] in contigs.keys():
            contigs[line[0]].assign_class_number(int(line[1]))
            count += 1
    print('> Processed {0}'.format(args.bins.name))

    print('> {0} unique bins found'.format(str(len(classes) - 1)))

    print('> {0} contigs assigned to a bin'.format(str(count)))

    # Color classes
    header_colors = ['%0 255\t255\t255']  # Default class is white
    rgb_tuples = rainbow_picker(len(classes) - 1)  # Skip first class
    for rgb_tuple in enumerate(rgb_tuples):
        color = '%{0} {1}\t{2}\t{3}'.format(classes[rgb_tuple[0] + 1],
                                            # Skip first class
                                            int(rgb_tuple[1][0]),
                                            int(rgb_tuple[1][1]),
                                            int(rgb_tuple[1][2]))
        header_colors.append(color)

    # Write .cls file
    print('> Writing {0}'.format(args.output.name))
    args.output.write('% {0}\n'.format(str(len(classes))))
    args.output.write('{0}\n'.format('\n'.join(header_colors)))
    chunk_ordered = {}
    for contig in contigs:
        for chunk in contigs[contig].chunk_numbers:
            chunk_ordered[chunk] = contigs[contig].class_number
    for chunk in sorted(chunk_ordered.keys()):
        args.output.write('{0}\t{1}\n'
                          .format(str(chunk),
                                  str(chunk_ordered[chunk])))
    args.output.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.
                                     RawDescriptionHelpFormatter)
    parser.add_argument('--bins', metavar='BINS file',
                        required=True,
                        type=x_reader,
                        help='file containing bins from anvi_converter')
    parser.add_argument('--names', metavar='NAMES file',
                        required=True,
                        type=x_reader,
                        help='NAMES file from ESOM')
    parser.add_argument('--output', metavar='OUT file',
                        required=True,
                        type=argparse.FileType('w'),
                        help='Output file to write, ".cls" will be added')
    args = parser.parse_args()

    main(args)

    sys.exit(0)
