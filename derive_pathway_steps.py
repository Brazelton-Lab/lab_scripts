#! /usr/bin/env python2

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
import os
import pycyc
import re
import stat
from subprocess import CalledProcessError, Popen
import sys
import time

__author__ = 'Alex Hyer'
__email__ = 'theonehyer@gmail.com'
__license__ = 'GPLv3'
__maintainer__ = 'Alex Hyer'
__status__ = 'Alpha'
__version__ = '0.0.1a12'


# This method is literally just the Python 3.5.1 which function from the
# shutil library in order to permit this functionality in Python 2.
# Minor changes to style were made to account for indentation.
def which(cmd, mode=os.F_OK | os.X_OK, path=None):
    """Given a command, mode, and a PATH string, return the path which
    conforms to the given mode on the PATH, or None if there is no such
    file.
    `mode` defaults to os.F_OK | os.X_OK. `path` defaults to the result
    of os.environ.get("PATH"), or can be overridden with a custom search
    path.
    """

    # Check that a given file can be accessed with the correct mode.
    # Additionally check that `file` is not a directory, as on Windows
    # directories pass the os.access check.
    def _access_check(fn, mode):
        return (os.path.exists(fn) and os.access(fn, mode)
                and not os.path.isdir(fn))

    # If we're given a path with a directory part, look it up directly
    # rather than referring to PATH directories. This includes checking
    # relative to the current directory, e.g. ./script
    if os.path.dirname(cmd):
        if _access_check(cmd, mode):
            return cmd
        return None

    if path is None:
        path = os.environ.get("PATH", os.defpath)
    if not path:
        return None
    path = path.split(os.pathsep)

    if sys.platform == "win32":
        # The current directory takes precedence on Windows.
        if not os.curdir in path:
            path.insert(0, os.curdir)
        # PATHEXT is necessary to check on Windows.
        pathext = os.environ.get("PATHEXT", "").split(os.pathsep)
        # See if the given file matches any of the expected path
        # extensions. This will allow us to short circuit when given
        # "python.exe". If it does match, only test that one, otherwise
        # we have to try others.
        if any(cmd.lower().endswith(ext.lower()) for ext in
               pathext):
            files = [cmd]
        else:
            files = [cmd + ext for ext in pathext]
    else:
        # On other platforms you don't have things like PATHEXT to tell you
        # what file suffixes are executable, so just pass on cmd as-is.
        files = [cmd]

    seen = set()
    for dir in path:
        normdir = os.path.normcase(dir)
        if not normdir in seen:
            seen.add(normdir)
            for thefile in files:
                name = os.path.join(dir, thefile)
                if _access_check(name, mode):
                    return name
    return None


def main(args):
    """Run program

        Args:
             args (NameSpace): ArgParse arguments controlling program flow
    """

    print('>>> Hi, I am DPS (Derive Pathway Steps).')
    print('>>> I will be analyzing pathways for you today.')
    print('>>> I am using the {0} database as per your command.'
          .format(args.database))

    if args.database == 'metacyc':

        # Obtain executable
        pathway_tools = which('pathway-tools')
        if pathway_tools is None:
            raise EnvironmentError('I cannot find pathway-tools: please '
                                   'specify -e')
        else:
            print('>>> I found pathway-tools: {0}'.format(pathway_tools))

        # Start pathway-tools daemon
        while True:
            print('>>> Starting Pathway Tools LISP Daemon.')
            pid = Popen([pathway_tools, '-lisp', '-api'],
                        stderr=open(os.devnull, 'w'),
                        stdout=open(os.devnull, 'w'))
            print('>>> Let\'s give it five seconds to spawn.')
            time.sleep(5)
            if os.path.exists('/tmp/ptools-socket') and \
                    stat.S_ISSOCK(os.stat('/tmp/ptools-socket').st_mode):
                print('>>> The daemon is is up!')
                break
            else:
                print('>>> The daemon took too long to boot. :(')
                print('>>> This makes me sad, so I will kill it.')
                pid.kill()
                print('>>> Let\'s wait five seconds for it to die!')
                time.sleep(5)
                pid.poll()
                if pid.returncode is None:
                    raise CalledProcessError('Pathway Tools won\'t die!')
                else:
                    print('>>> The daemon is dead!')
                    print('>>> I miss it. :( I\'m going to try again. :)')

        # Connect to daemon
        try:
            metacyc = pycyc.open('meta')
        except IOError:
            print('>>> I cannot connect to Pathway Tools Daemon.')
            print('>>> Here is the original error message:')
            raise
        else:
            print('>>> I have connected to the Pathway Tools Daemon.')
            print('>>> Phenomenal cosmic powers! Itty bitty memory footprint!')

        # Index genes file
        print('>>> Indexing {0}.'.format(args.reactions_file.name))
        reactions_to_genes = {}
        start_time = time.time()
        for line in args.genes_file:
            parts = line.strip().split()
            reactions_to_genes[parts[0]] = (parts[1], [' '.join(parts[2:])])
        end_time = time.time()
        print('>>> I indexed {0} reactions in {1} seconds.'
              .format(str(len(reactions_to_genes)),
                      str(end_time - start_time)))
        print('>>> I\'m so fast.')

        # Index UniRef mapping file
        print('>>> Indexing {0}.'.format(args.uniref_file.name))
        unirefs = {}
        start_time = time.time()
        for line in args.uniref_file:
            parts = line.strip().split()
            unirefs[parts[0]] = (parts[1], parts[3])
        end_time = time.time()
        print('>>> I indexed {0} reactions in {1} seconds.'
              .format(str(len(unirefs)),
                      str(end_time - start_time)))
        print('>>> I\'m quick!')

        # Index all pathways by name
        print('>>> Time to index all the pathways from Metacyc.')
        pathways = {}
        start_time = time.time()
        for frame in metacyc.all_pathways():
            pathways[frame.common_name] = frame
        end_time = time.time()
        print('>>> I indexed {0} pathways in {1} seconds.'
              .format(str(len(pathways)), str(end_time - start_time)))
        print('>>> Aren\'t you proud of me?')

        # Obtain pathway of interest
        print('>>> Time to do some science!')
        print('>>> Note: you can input all or part of a pathway name.')
        print('>>> Type "q" for input at any time to exit the program.')

        while True:  # Rest of program runs in a loop until user ends it
            possiblities = {}
            user_input = raw_input('>>> Enter a pathway: ')
            if user_input.lower() == 'q':
                print('>>> Shutdown sequence initiated.')
                break
            for name, frame in pathways.items():
                if user_input in name:
                    possiblities[name] = frame
            if len(possiblities) == 0:
                print('>>> I couldn\'t find any pathways matching your '
                      'request.')
                print('>>> Try ')
            print('>>> I found {0} pathways matching your request.'
                  .format(str(len(possiblities))))

            shutdown = False
            pathway = None
            while True:
                print('>>> Here are possible pathways:')
                max_entry = len(possiblities) - 1
                for possibility in enumerate(possiblities):
                    print('{0}: {1}'.format(str(possibility[0]),
                                            possibility[1].common_name))
                path_num = raw_input('>>> Select a pathway: ')
                if path_num.lower() == 'q':
                    shutdown = True
                    break
                else:
                    try:
                        path_num = int(path_num)
                    except ValueError:
                        print('>>> Your answer is not an integer.')
                        print('>>> I only understand integers.')
                        print('>>> Please correct.')
                        continue
                    if path_num > max_entry or path_num < 0:
                        print('>>> {0} is not a valid pathway.')
                        print('>>> Valid pathways are: {0}.'.format(' '.join(
                                [str(i) for i in range(max_entry)])))
                        print('>>> Try again.')
                        continue
                    pathway = possiblities[possiblities.keys[path_num]]
                    print('>>> You selected: {0}.')
                    print('>>> Neat! I\'ll analyze it now.')
                    break
            if shutdown is True:
                print('>>> Shutdown sequence initiated.')
                break

            # Add genes and abundances to pathway reactions
            print('>>> Collecting reactions in pathway.')
            rxns = [str(rxn) for rxn in pathway.reaction_list]
            print('>>> Analyzing pathway for key genes.')
            key_rxns = [str(key) for key in pathway.key_reactions]

            print('>>> Acquiring gene families for each reaction from {0}. '
                  .format(args.reactions_file.name))
            reactions = {}
            for rxn in rxns:
                if rxn in reactions_to_genes.keys():
                    name, uniref_list = reactions_to_genes[rxn]
                    if rxn in key_rxns:
                        name += '*'
                    name = name + ' (' + uniref_list[0] + ')'
                    reactions[name] = {}
                    for uniref in uniref_list[1].split():
                        reactions[name][uniref] = {}

            print('>>> Adding abundances from {0}.'
                  .format(args.abundance_file.name))
            for line in args.abundance_file:
                gene, abundance = line.strip().split('\t')
                if gene in unirefs:
                    uniref, name = unirefs[gene]
                    for rxn in reactions.keys():
                        if uniref in reactions[rxn].keys():
                            reactions[rxn][uniref] = (abundance, name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.
                                     RawDescriptionHelpFormatter)

    subparsers = parser.add_subparsers(title='Database',
                                       dest='database')

    metacyc = subparsers.add_parser('metacyc',
                                    help='Analyze MetaCyc Database')
    metacyc.add_argument('abundance_file',
                         metavar='Abundance File',
                         type=argparse.FileType('r'),
                         help='TSV containing gene ID and abundance columns')
    metacyc.add_argument('reactions_file',
                         metavar='Reactions File',
                         type=argparse.FileType('r'),
                         help='metacyc1 file mapping Unirefs to reactions')
    metacyc.add_argument('uniref_file',
                         metavar='UniRef ID File',
                         type=argparse.FileType('r'),
                         help='ID Mapping file mapping UniRefs to genes')
    metacyc.add_argument('-e', '--executable',
                         default=None,
                         type=str,
                         help='pathways-tree executable if not in PATH')
    args = parser.parse_args()

    main(args)

    sys.exit(0)
