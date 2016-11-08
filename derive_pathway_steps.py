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
import os
import stat
from subprocess import CalledProcessError, Popen
import sys
import time

__author__ = 'Alex Hyer'
__email__ = 'theonehyer@gmail.com'
__license__ = 'GPLv3'
__maintainer__ = 'Alex Hyer'
__status__ = 'Alpha'
__version__ = '0.0.1a11'


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

        # Start pathway-tools server
        while True:
            print('>>> Starting Pathway Tools LISP Daemon.')
            pid = Popen([pathway_tools, '-lisp', '-api'],
                        stderr=open(os.devnull), stdout=open(os.devnull))
            print('>>> Let\'s give it five seconds.')
            time.sleep(5)
            if os.path.exists('/tmp/ptools-socket') and \
                    stat.S_ISSOCK(os.stat('/tmp/ptools-socket').st_mode):
                print('>>> The daemon is is up!')
                break
            else:
                print('>>> The daemon took too long to boot. :(')
                print('>>> This makes me sad so I will kill it.')
                pid.kill()
                print('>>> Let\'s wait five seconds for it to die!')
                time.sleep(5)
                pid.poll()
                if pid.returncode is None:
                    raise CalledProcessError('Pathway Tools won\'t die!')
                else:
                    print('>>> The daemon is dead!')
                    print('>>> I miss it. :( I\'m going to try again. :)')


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

