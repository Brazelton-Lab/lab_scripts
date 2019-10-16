#! /usr/bin/env python

"""
given a .gff file and a set of fasta files containing contigs represented in the gff file, divide the gff file into separate gff files, one for each of the provided contig fasta files to be used to extract separate gff files for bins that are subsets of a larger assembly must be run from directory containing fasta files.

usage:
python get-gff.py assembly.gff fasta

Copyright:

    get-gff  divides a gff file into separate gff files, one for each of the provided contig fasta files

    Copyright (C) 2019  William Brazelton

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

# define input files
import sys
gff_file = sys.argv[1]
extension = '*' + sys.argv[2]

# check that files are in proper order in command line
if gff_file[-4:] == '.gff': pass
else:
	print 'Please enter files in this order: assembly.gff extension-of-fasta-files-with-contigs'
	sys.exit()


# create dictionary of fasta filename = contig names in that fasta file
c = {}
l = {}
import glob
for filename in glob.glob(extension):
	with open(filename) as contigs:
		count = 0
		for line in contigs: 
			if line[0] == '>': 
				count = count + 1
				header = line.strip('>')
				contig_id = header.strip('\n')
				if filename in c: c[filename].append(contig_id)
				else: c[filename] = list(tuple([contig_id]))
		l[filename] = count

# write each component of the gff file to the appropriate file			
for filename in c:
	newfilename = filename.replace(sys.argv[2],'gff')
	with open(newfilename,'w') as newfile:
		newfile.write('##gff-version 3\n')
		with open(gff_file) as gff:
			scount = 0
			gcount = 0
			status = 'no' 			# this switch controls whether the "# Model Data" header line is written
			for line in gff:
				if "# Model" in line:
					if status == 'yes':
						newfile.write(line)
						status = 'no'
				if "# Sequence" in line:
					cols = line.split(';')
					seq = cols[2].split('seqhdr=')
					seqname = seq[1].replace('"','')
					seqname = seqname.strip('\n')
					if seqname in c[filename]:
						newfile.write(line)
						scount = scount + 1
						status = 'yes'
					else: status = 'no'
				else:  
					cols = line.split('\t')
					contig_name = cols[0]
					if contig_name in c[filename]: 
						newfile.write(line)
						gcount = gcount + 1
						status = 'no'
	numfasta = l[filename]
	print str(numfasta) + ' FASTA entries in ' + filename
	print str(scount) + ' sequence header lines written to ' + newfilename
	print str(gcount) + ' gff3 format lines written to ' + newfilename
