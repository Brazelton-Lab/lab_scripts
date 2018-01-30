#! /usr/bin/env python

"""
edits mothur taxonomy file
Transfers last name that is not "unclassified" or "uncultured" to "unclassified" or "uncultured" assignment. Also removes numbers in parentheses.

Copyright:

    taxonomy_edit  Cleans up a mothur taxonomy file

    Copyright (C) 2016  William Brazelton

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.﻿
"""

import sys
infilename = sys.argv[1]
outfilename = infilename + '.renamed.txt'

with open(infilename) as infile:
        for line in infile:
                columns = line.split('\t')
                with open(outfilename,'a') as outfile: outfile.write(columns[0])

                tax = columns[1]
                if "unclassified" in tax:
                        names = tax.split(';')
                        for name in names:
                                if name == "unclassified":
                                        newname = name + '_' + prevgoodname
                                        with open(outfilename,'a') as outfile: outfile.write('\t' + newname)
                                elif name == "uncultured":	# because both unclassified and uncultured could be in same taxonomy line
                                        newname = name + '_' + prevgoodname
                                        with open(outfilename,'a') as outfile: outfile.write('\t' + newname)
                                else:
                                     	name = name.split('(')
                                        name = name[0]
                                        if name == "uncultured":        # because some "uncultured"s have numbers and parentheses and are not caught above
                                                newname = name + '_' + prevgoodname
                                                with open(outfilename,'a') as outfile: outfile.write('\t' + newname)
                                        else:
                                             	with open(outfilename,'a') as outfile: outfile.write('\t' + name)
                                                prevgoodname = name

                elif "uncultured" in tax:
                        names = tax.split(';')
                        for name in names:
                                if name == "uncultured":
                                        newname = name + '_' + prevgoodname
                                        with open(outfilename,'a') as outfile: outfile.write('\t' + newname)
                                else:
                                     	name = name.split('(')
                                        name = name[0]
                                        if name == "uncultured":        # because some "uncultured"s have numbers and parentheses and are not caught above
                                                newname = name + '_' + prevgoodname
                                                with open(outfilename,'a') as outfile: outfile.write('\t' + newname)
                                        else:
                                             	with open(outfilename,'a') as outfile: outfile.write('\t' + name)
                                                prevgoodname = name
                else:
                     	names = tax.split(';')
                        for name in names:
                                name = name.split('(')
                                name = name[0]
                                with open(outfilename,'a') as outfile: outfile.write('\t' + name)
