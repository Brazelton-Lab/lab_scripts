#! /usr/bin/env python
# counts number of "mapped read pairs" = "mapped fragments"
# this is NOT the same thing as number of "mapped paired reads" = "mapped reads"
# uses the bowtie2 log file where "mapped read pairs" = "were paired" - "pairs aligned 0 times concordantly or discordantly"
# i.e. the number of read pairs subtracted by the number of read pairs that did NOT map (concordantly or discordantly)
# usage:
# python count-mapped-fragments.py
# will print to screen the name of every file ending in ".log" and the number calculated from the equation described above

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
