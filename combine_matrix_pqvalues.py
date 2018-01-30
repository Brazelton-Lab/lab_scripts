#! /usr/bin/env python


"""
extract pairwise distances from matrix and combine with pairs in columns from another file
use with cor.mat file from rcor.test function in R and output file from QVALUE

usage: python combine_matrix_pqvalues.py matrix_filename qvalue_filename


Copyright:

    combine_matrix_pqvalues.py Extract pairwise distances and combine w/ rcor

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

import math
import sys
file1name = sys.argv[1]
file2name = sys.argv[2]

file1 = open(file1name)
file2 = open(file2name)
outfile = open('your_cors_and_pqvalues.txt', 'a')
outfile.write('var1\tvar2\tcor_coeff\ttype\tp-value\tq-value\n')

# collect names and rows from matrix
matrix = file1.read()
rows = matrix.split('\n')
name_list = []
row_list = []
for name in (rows[0].split(','))[1:]:
	name_list.append(name)
for line in rows[1:-1]:
	columns = line.split(',')
	row = []
	for column in columns[1:]:
		row.append(column)			# creates list of items in that row except the row header
	row_list.append(row)			# creates tuple of lists of entries

# iterate through p-values and q-values and write the appropriate value from the stored tuple of rows
row_count = 0
value_count = 1
num_rows = len(name_list)
for pair in file2.xreadlines():
	if pair[0] == 'p': pass
	elif len(pair) < 2: pass
	else:
		pqvalues = pair.split(' ')
		outfile.write(name_list[row_count])
		outfile.write('\t')
		outfile.write(name_list[value_count])
		outfile.write('\t')

		current_row = row_list[row_count]
		current_entry = current_row[value_count]
		if float(current_entry) > 0: sign = 'pos'
		elif float(current_entry) < 0: sign = 'neg'
		elif float(current_entry) == 0: sign = 'neutral'
		current_entry = abs(float(current_entry))
		outfile.write(str(current_entry))
		outfile.write('\t')

		try: outfile.write(sign)
		except: outfile.write('error')
		outfile.write('\t')

		outfile.write(str(pqvalues[0]))
		outfile.write('\t')
		outfile.write(str(pqvalues[1]))
		value_count = value_count + 1
		if value_count == num_rows:
			row_count = row_count + 1
			value_count = row_count + 1
		if row_count + 1 == len(name_list): break

file1.close()
file2.close()
outfile.close()
