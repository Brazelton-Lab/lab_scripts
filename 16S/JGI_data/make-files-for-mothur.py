#! /usr/bin/env python
# from folder of .fastq files, will generate one file listing the sample name, file1name, file2name
# assumes filenames have structure: name.1.fastq and name.2.fastq

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