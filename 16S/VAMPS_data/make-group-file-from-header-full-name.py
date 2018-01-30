#! /usr/bin/env python

"""
sort sequences in standard.fasta according to bar code
files saved with filenames according to barcodes_distribution.txt
this version has the name parsing commented out in order to include the full sample name

Copyright:

    make-group-file-from-header-full-name  sort sequences in standard.fasta according to bar code

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
infilename = sys.argv[1]

outfilename = infilename + '.group'
outfile = open(outfilename,'a')
                
fastas = open(infilename)
for line in fastas:
        #print line
        if line[0] == '>':
                header = line.strip('\n')
                header = header.split(' ')
                header = header[0]
                words = line.split('|')
                sample = words[2]
                
                #fasta_sample = words[2].split('_')
                #if len(fasta_sample) < 3: sample = fasta_sample
                #else: sample = fasta_sample[1:-1]

                outfile.write(header.replace('>',''))
                outfile.write('\t')
                outfile.write(sample)
                
                #outfile.write(sample[0])
                #for i in sample[1:]:
                #        outfile.write('_')
                #        outfile.write(i)
                outfile.write('\n')
                
        else: pass
fastas.close()  
outfile.close()
