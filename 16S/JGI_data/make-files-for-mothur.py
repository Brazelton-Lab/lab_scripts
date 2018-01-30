#! /usr/bin/env python


"""
from folder of .fastq files, will generate one file listing the sample name, file1name, file2name
assumes filenames have structure: name.1.fastq and name.2.fastq

Copyright:

    make-files-for-mothur  from folder of .fastq files, will generate one file listing the sample name, file1name, file2name

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

import os
import glob

outfile=open('JGI.files','a')

previous=''
path = r'./'
for filename in glob.glob(os.path.join(path, '*.fastq')):
	filename = filename.replace('./','')
	print filename
	split = filename.split('.')
	if split[0] == previous: pass
	else:
		previous = split[0]
		outfile.write(split[0])
		outfile.write('\t')
		outfile.write(filename)
		outfile.write('\t')
		outfile.write(filename.replace('1.fastq','2.fastq'))
		outfile.write('\n')
		
outfile.close()		
