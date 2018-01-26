#! /usr/bin/env python

"""
Copyright:

    MSU-16S-arch-usearch8-wrapper  Reformats a file's filename and calls usearch8 on it.

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
import re
import subprocess
import sys

dir = sys.argv[1]
if dir == '.':
    dir = os.getcwd()
elif dir.endswith('\/'):
    dir = dir[:-1]

for file in glob.glob(dir + os.sep + '*'):
    if '_R2_' in file:
        pass
    else:
        rfile = file.replace('_R1_', '_R2_')
        fastqOut = re.sub('_R1_.*', '.merged.fastq', file)
        firstLog = fastqOut.replace('.fastq', '.log')
        fastaOut = fastqOut.replace('.fastq', '.filtered.fasta')
        secondLog = firstLog.replace('.log', '.filtered.log')
        subprocess.call(['usearch8', '-fastq_mergepairs', file,
                         '-reverse', rfile, '-fastqout', fastqOut,
                         '-fastq_minovlen', '30', '--log', firstLog])
        subprocess.call(['usearch8', '-fastq_filter', fastqOut,
                         '-fastaout', fastaOut, '-fastq_maxee', '1',
                         '--log', secondLog])
