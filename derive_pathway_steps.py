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
__version__ = '0.0.1a14'


def print_list(lst, level=0):
    yield('    ' * (level - 1) + '+---' * (level > 0) + str(lst[0]))
    for l in lst[1:]:
        if type(l) is list:
            for i in print_list(l, level + 1):
                yield i
        else:
            yield('    ' * level + '+---' + str(l))


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
                                   'specify -e.')
        else:
            print('>>> I found pathway-tools: {0}.'.format(pathway_tools))

        # Start pathway-tools daemon
        while True:
            print('>>> Summoning Pathway Tools LISP Daemon.')
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
        for line in args.reactions_file:
            parts = line.strip().split()
            reactions_to_genes[parts[0]] = (parts[1], parts[2:])
        end_time = time.time()
        print('>>> I indexed {0} reactions in {1} seconds.'
              .format(str(len(reactions_to_genes)),
                      str(end_time - start_time)))
        print('>>> I\'m so fast.')

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
            possibilities = {}
            user_input = raw_input('>>> Enter a pathway: ')
            if user_input.lower() == 'q':
                break
            for name, frame in pathways.items():
                if user_input in name:
                    possibilities[name] = frame
            if len(possibilities) == 0:
                print('>>> I couldn\'t find any pathways matching your '
                      'request.')
                print('>>> Try ')
            print('>>> I found {0} pathways matching your request.'
                  .format(str(len(possibilities))))

            shutdown = False
            pathway = None
            while True:
                print('>>> Here are possible pathways:')
                max_entry = len(possibilities) - 1
                for possibility in enumerate(possibilities.items()):
                    print('{0}: {1}'.format(str(possibility[0]),
                                            possibility[1][1].common_name))
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
                        print('>>> {0} is not a valid pathway.'
                              .format(str(path_num)))
                        print('>>> Valid pathways are: {0}.'.format(' '.join(
                                [str(i) for i in range(max_entry + 1)])))
                        print('>>> Try again.')
                        continue
                    pathway = possibilities[possibilities.keys()[path_num]]
                    print('>>> You selected: {0}.'.format(pathway.common_name))
                    print('>>> Neat! I\'ll analyze it now.')
                    break
            if shutdown is True:
                break

            # Add genes and abundances to pathway reactions
            print('>>> Collecting reactions in pathway.')
            rxns = [str(rxn) for rxn in pathway.reaction_list]
            print('>>> Analyzing pathway for key reactions.')
            if pathway.key_reactions is not None:
                key_rxns = [str(key) for key in pathway.key_reactions]
                for rxn in enumerate(rxns):
                    if rxn[1] in key_rxns:
                        rxns[rxn[0]] = rxn[1] + '*'

            print('>>> Acquiring gene families for each reaction from {0}.'
                  .format(args.reactions_file.name))
            reactions = {}
            for rxn in rxns:
                rxn_name = re.sub('\*$', '$', rxn)
                if rxn_name in reactions_to_genes.keys():
                    ec, uniref_list = reactions_to_genes[rxn_name]
                    rxn_name = rxn_name + ' (' + ec + ')'
                    reactions[rxn_name] = {}
                    for uniref in uniref_list:
                        reactions[rxn_name][uniref] = 0.0

            print('>>> Adding abundances from {0}.'
                  .format(args.abundance_file.name))
            for line in args.abundance_file:
                uniref, abundance = line.strip().split('\t')
                if float(abundance) > 0.0:
                    for rxn in reactions.keys():
                        if uniref in reactions[rxn].keys():
                            reactions[rxn][uniref] = abundance

            print('>>> Removing unused gene families.')
            for rxn in reactions.keys():
                for uniref in reactions[rxn].keys():
                    if reactions[rxn][uniref] == 0.0:
                        del reactions[rxn][uniref]
            for rxn in reactions.keys():
                if reactions[rxn] == {}:
                    reactions[rxn] = 'None\tN/A'
                    continue

            rxn_list = [pathway.common_name]
            for rxn in reactions.keys():
                if reactions[rxn] == 'None\tN/A':
                    temp = [rxn, ['None\tN/A']]
                    rxn_list.append(temp)
                elif type(reactions[rxn]) is dict:
                    temp = [rxn]
                    sub_temp = []
                    for uniref in reactions[rxn].keys():
                        sub_temp.append('{0}\t{1}'.format(uniref,
                                        str(reactions[rxn][uniref])))
                    temp.append(sub_temp)
                    rxn_list.append(temp)

            # Print output
            print('>>> I\'ve finished analyzing everything!')
            print('>>> Here it is (asterisks represent key reactions):')
            rxn_print = [rxn for rxn in print_list(rxn_list)]
            for rxn in rxn_print:
                print(rxn)

            # Save output
            print('>>> What file would you like me to save this to?')
            print('>>> Type "n" if you don\'t want to save this output.')
            while True:
                out_file = raw_input('>>> File: ')
                if out_file.lower() != 'n':
                    try:
                        with open(out_file, 'w') as out_handle:
                            for rxn in rxn_print:
                                out_handle.write(rxn + os.linesep)
                        print('>>> Output written to {0}.'.format(out_file))
                        break
                    except IOError as error:
                        print('>>> I could not write to {0}.'.format(out_file))
                        print('>>> Original error:')
                        print(error)
                        print('>>> Let\'s try again (enter "n" to skip).')
                else:
                    print(dir(pathway))
                    break

    # Shutdown
    print('>>> Shutdown sequence initiated.')


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
    metacyc.add_argument('-e', '--executable',
                         default=None,
                         type=str,
                         help='pathways-tree executable if not in PATH')
    args = parser.parse_args()

    main(args)

    sys.exit(0)
