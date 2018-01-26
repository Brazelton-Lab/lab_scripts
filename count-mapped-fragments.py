#! /usr/bin/env python


"""
counts number of "mapped read pairs" = "mapped fragments"
this is NOT the same thing as number of "mapped paired reads" = "mapped reads"
uses the bowtie2 log file where "mapped read pairs" = "were paired" - "pairs aligned 0 times concordantly or discordantly"
i.e. the number of read pairs subtracted by the number of read pairs that did NOT map (concordantly or discordantly)

usage:
python count-mapped-fragments.py
will print to screen the name of every file ending in ".log" and the number calculated from the equation described above

Copyright:

    count-mapped-fragments.py Add up mapped fragments from Bowtie2 output
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
for filename in glob.glob('*.log'):
	with open(filename) as f:
		pairs = None
		not_mapped = None
		for line in f:
			if "were paired; of these:" in line:
				line = line.lstrip()
				line = line.split()
				pairs = int(line[0])
			elif "pairs aligned 0 times concordantly or discordantly; of these:" in line:
				line = line.lstrip()
				line = line.split()
				not_mapped = int(line[0])
			else: pass
		try:
			mapped_pairs = pairs - not_mapped
			print mapped_pairs,
			print filename
		except: pass
