#! /usr/bin/env python

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
                         '-fastq_minovlen', '60', '--log', firstLog])
        subprocess.call(['usearch8', '-fastq_filter', fastqOut,
                         '-fastaout', fastaOut, '-fastq_maxee', '1',
                         '--log', secondLog])
