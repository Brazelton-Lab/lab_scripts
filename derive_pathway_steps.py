#! /usr/bin/env python

"""

Copyright:

    derive_pathway+steps.py Obtain gene list from pathway databases
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
from time import time
import sys

__author__ = 'Alex Hyer'
__email__ = 'theonehyer@gmail.com'
__license__ = 'GPLv3'
__maintainer__ = 'Alex Hyer'
__status__ = 'Alpha'
__version__ = '0.0.1a4'


def main(args):
    """Run program

        Args:
             args (NameSpace): ArgParse arguments controlling program flow
    """

    print('>>> Hi, I\'m DPS (Derive Pathway Steps)')
    print('>>> I will be analyzing pathways for you today')
    print('>>> I am using the {0} database as per your command'
          .format(args.database))

    # Index genes file
    print('>>> I am indexing {0}'.format(args.genes_file.name))
    print('>>> I\'ll try to be super fast')
    rxn_index = {}
    index = 0
    start_time = time()
    while True:
        line = args.genes_file.readline()
        rxn = line.split('\t').strip()[0]
        rxn_index[rxn] = index
        index = args.genes_file.tell()
    end_time = time()
    print('>>> I indexed {0} reactions in {1} seconds'
          .format(str(len(rxn_index.keys())), str(end_time - start_time)))
    print('>>> See how quickly I did that?')

    # Find possible pathways
    pathway_index = {}
    index = 0
    start_time = time()
    while True:
        line = args.pathways_file.readline()
        pathway = line.split('\t').strip()[0]
        if args.pathway in pathway:
            pathway_index[pathway] = index
        index = args.pathways_file.tell()
    end_time = time()
    print('>>> I indexed {0} possible pathways in {1} seconds'
          .format(str(len(pathway_index.keys())), str(end_time - start_time)))
    print('>>> I hope they\'re relevant to your bioinformatic desires')
    print('>>> Select the pathway you want me to analyze')
    answer = -1
    index = None
    while True:
        print('>>> Possible Pathways (select a number)')
        for item in enumerate(pathway_index.items()):
            print('{0}: {1}'.format(str(item[0]), item[1][0]))
        answer = input('>>> Tell me the what pathway to analyze: ')
        try:
            index = [i[1] for i in pathway_index.items()][answer]
            break
        except IndexError:
            print('>>> I don\'t know what pathway {0} is. Please clarify.')
    sys.exit(0)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.
                                     RawDescriptionHelpFormatter)
    parser.add_argument('output_file',
                        metavar='Output File',
                        type=argparse.FileType('w'),
                        help='output file for results')

    subparsers = parser.add_subparsers(title='Database',
                                       dest='database')

    metacyc = subparsers.add_parser('metacyc',
                                    metavar='MetaCyc',
                                    help='Analyze MetaCyc Database')
    metacyc.add_argument('pathways_file',
                         metavar='Pathways File',
                         type=argparse.FileType('r'),
                         help='metacyc2 file: connects enzymatic reactions to '
                              'pathways')
    metacyc.add_argument('genes_file',
                         metavar='Genes File',
                         type=argparse.FileType('r'),
                         help='metacyc1 file: connects genes to enzymatic '
                              'reactions')
    metacyc.add_argument('uniref_file',
                         metavar='UniRef ID File',
                         type=argparse.FileType('r'),
                         help='ID Mapping file mapping UniRefs to genes')
    metacyc.add_argument('pathway',
                         metavar='Pathway',
                         type=str,
                         help='Name of pathway to analyze (can be partial)')
    args = parser.parse_args()

    main(args)

    sys.exit(0)

