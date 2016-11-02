#! /usr/bin/env python

# report column with the highest number
# assumes comma-separated-values and unix linebreaks
# usage:
# python find-max-column.py final.an.0.03.subsample.shared.nosingles.transposed.csv

"""
Copyright:

    find-max-column.py Report column in CSV with highest number
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
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import sys
filename = sys.argv[1]
file = open(filename)
outfilename = filename.replace('.csv','.maxcolumn.csv')
outfile = open(outfilename, 'a')

for line in file:
        line = line.strip('\n')
        headers = line.split(',')
        headers = headers[1:]
        break

for line in file:
        line = line.strip('\n')
        columns = line.split(',')
        OTU = columns[0]
        numbers = []
        for i in columns[1:]: numbers.append(float(i))
        maxcolumn = numbers.index(max(numbers))
        outfile.write(OTU)
        outfile.write(',')
        outfile.write(headers[maxcolumn])
        outfile.write(',')
        outfile.write(str(max(numbers)))
        outfile.write('\n')

file.close()
outfile.close()
