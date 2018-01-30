#! /usr/bin/env python

"""
report length of each FASTA sequence and report average

Copyright:

    FASTA_size_and_number_loop  report length of each FASTA sequence and report average

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

import glob
import os
path = r'./'
for fileName in glob.glob(os.path.join(path, '*.fna')):
	infile = open(fileName)
	length_list = []
	length = 0
	for line in infile:
		if line[0] == '>': 
			if length == 0: pass
			else:
				length_list.append(length)
				length = 0
		else:
			length = length + len(line)
	length_list.append(length)		
				
	total = 0
	for entry in length_list:
		total = int(total) + int(entry)
	avg = float(total) / float(len(length_list))	
	print fileName
	print len(length_list),
	print " FASTA sequences"
	print total,
	print " total bps"
	print avg,
	print " average length"
	print min(length_list),
	print " is shortest sequence"
	print max(length_list),
	print " is longest sequence"


	
