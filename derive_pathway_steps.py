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
import os
import re
import ruamel.yaml
import sys

__author__ = 'Alex Hyer'
__email__ = 'theonehyer@gmail.com'
__license__ = 'GPLv3'
__maintainer__ = 'Alex Hyer'
__status__ = 'Alpha'
__version__ = '0.0.1a11'


class LispParserError(Exception):
    pass

token_pattern = re.compile(r'\s*([()]|"[^"]*"|[^"\s()]*)(.*)', re.DOTALL)

class LispParser():

    def parse_string(self,expression):
        """
        Tokenize and parse a string as a LISP expression.
        """
        try:
            return self.parse_token_generator(self.token_generator(expression))
        except StopIteration:
            # This should occur only if the expression is malformed
            raise LispParserError('Malformed or unsupported LISP expression.')

    def token_generator(self,expression):
        """ etc """

        while expression:
            token, expression = token_pattern.match(expression).groups()
            if token:
                yield token
            else:
                # Handle pathological cases like expression='"', which
                # gives token = '' and expression = '"' again
                if expression:
                    raise LispParserError('Failure to parse LISP string')
                break

    def parse_token_generator(self, token_generator):
        """
        Parse a generator of LISP token strings.
        """
        token = next(token_generator)

        if token == '(':
            elements = [token]
            while elements[-1] != ')':
                elements.append(self.parse_token_generator(token_generator))
            return elements[1:-1]
        elif token == ')':
            return ')'
        else:
            return self.parse_single_token(token)

    def parse_single_token(self, token):
        """
        Parse a single LISP token string to a number, boolean, None, or string.
        """

        if token.startswith('"'):
            return token.strip('"')

        if token == 'T':
            return True
        elif token == 'NIL':
            return None

        try:
            return int(token)
        except ValueError:
            pass
        try:
            return float(token)
        except ValueError:
            pass

        return token


class Pathway(object):
    """A class to store the structure of a MetaCyc tree

    Attributes:
        name (str): name of pathway

        child_nodes (dict): dictionary of Reaction consisting of all nodes one
                            level down from current node, key is child node's
                            name and value is pointer to node
    """

    def __init__(self, name):
        """Initialize attributes to store pathway structure"""

        self.name = name
        self.child_nodes = {}

    def add_child_node(self, child_name):
        """Add child node to child_nodes with self as parent

        Args:
            child_name (str): Reaction name of child node
        """

        self.child_nodes[child_name] = Reaction(child_name, self)

    def get_tree_structure(self):
        """Return lists of all possible reaction paths in pathway as objects

        Returns:
             list: list of lists of all possible branches as objects
        """

        branches = []
        for child in self.child_nodes.values():
            for branch in child.gather_children():
                branches.append(branch)

        return branches

    def print_tree_structure(self):
        """Return lists of all possible reaction paths in pathway as strings

        Returns:
             list: list of lists of all possible branches as strings
        """

        branches = []
        for child in self.child_nodes.values():
            for branch in child.gather_children_names():
                branches.append(branch)

        return branches


class Reaction(Pathway):
    """A class to act as a node on a MetaCyc pathway tree

    Attributes:
        name (str): Name of node

        parent(Pathway|Reaction): Pathway or Reaction node above current node
    """

    def __init__(self, name, parent):
        """Initialize node and call parent __init__"""

        super(Reaction, self).__init__(name)
        self.parent_name = parent.name
        self.parent_node = parent

    def gather_children(self):
        """Generator of lists containing each path branching off current node

        Yields:
            list: list of all children objects along a given branch
        """

        # Only return own name and cease iteration if node is terminal leaf
        if not self.child_nodes:
            yield [self]
        else:
            for child in self.child_nodes.values():
                branch = [self]
                for piece in child.gather_children():  # Recurse
                    yield branch + piece

    def gather_children_names(self):
        """Generator of lists containing each path branching off current node

        Yields:
            list: list of all children names along a given branch
        """

        # Only return own name and cease iteration if node is terminal leaf
        if not self.child_nodes:
            yield [self.name]
        else:
            for child in self.child_nodes.values():
                branch = [self.name]
                for piece in child.gather_children_names():  # Recurse
                    yield branch + piece


def metacyc_tree(raw_data, name):
    """Analyzes MetaCyc structured trees and parses into a tree

    Args:
        raw_data (str): raw structured format of tree

        name (str): name of pathway being analyzed

    Returns:
        Pathway: Pathway class modeling pathway
    """

    # TODO: Download PathwayTools
    pass


def main(args):
    """Run program

        Args:
             args (NameSpace): ArgParse arguments controlling program flow
    """

    print('>>> Hi, I\'m DPS (Derive Pathway Steps)')
    print('>>> I will be analyzing pathways for you today')
    print('>>> I am using the {0} database as per your command'
          .format(args.database))

    # Index genes database
    rxn_index = {}
    print('>>> I am indexing {0}'.format(args.genes_file.name))
    print('>>> I\'ll try to be super fast')
    print('>>> Buckle up!')
    index = 0
    start_time = time()
    line = args.genes_file.readline()
    while line:
        rxn = line.strip().split('\t')[0]
        rxn_index[rxn] = index
        index = args.genes_file.tell()
        line = args.genes_file.readline()
    end_time = time()
    print('>>> I indexed {0} reactions in {1} seconds'
          .format(str(len(rxn_index.keys())), str(end_time - start_time)))
    print('>>> See how quickly I did that?')
    print('>>> I\'ll attempt to save an index file to save computation '
          'next time')

    # Find possible pathways and index them
    print('>>> Finding pathways with {0} in {1}'
          .format(args.pathway, args.pathways_file.name))
    pathway_index = {}
    index = 0
    start_time = time()
    line = args.pathways_file.readline()
    while line:
        pathway = line.strip().split('\t')[0]
        if args.pathway in pathway:
            pathway_index[pathway] = index
        index = args.pathways_file.tell()
        line = args.pathways_file.readline()
    end_time = time()
    print('>>> I indexed {0} possible pathways in {1} seconds'
          .format(str(len(pathway_index.keys())), str(end_time - start_time)))
    print('>>> I hope they\'re relevant to your bioinformatic desires')
    print('>>> Select the pathway you want me to analyze')
    answer = -1
    index = None
    pathway_name = None
    while True:
        print('>>> Possible Pathways (select a number)')
        for item in enumerate(pathway_index.items()):
            print('{0}: {1}'.format(str(item[0]), item[1][0]))
        answer = input('>>> Tell me the what pathway to analyze: ')
        try:
            index = [i[1] for i in pathway_index.items()][answer]
            pathway_name = [i[0] for i in pathway_index.items()][answer]
            break
        except IndexError:
            print('>>> I don\'t know what pathway {0} is. Please clarify.'
                  .format(str(answer)))

    # Analyze structure of pathway
    print('>>> I will now analyze {0}'.format(pathway_name))
    args.pathways_file.seek(index)
    pathway_structure_raw = args.pathways_file.strip().split('\t')[1]
    pathway_tree = metacyc_tree(pathway_structure_raw, pathway_name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.
                                     RawDescriptionHelpFormatter)

    subparsers = parser.add_subparsers(title='Database',
                                       dest='database')

    metacyc = subparsers.add_parser('metacyc',
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
    metacyc.add_argument('output_file',
                         metavar='Output File',
                         type=argparse.FileType('w'),
                         help='output file for results')
    metacyc.add_argument('pathway',
                         metavar='Pathway',
                         type=str,
                         help='Name of pathway to analyze (can be partial)')
    args = parser.parse_args()

    main(args)

    sys.exit(0)

