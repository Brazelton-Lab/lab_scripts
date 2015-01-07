#! /usr/bin/env python
# runs prinseq-lite on all .fastq files in a folder

import os
import glob

path = r'./'
for filename in glob.glob(os.path.join(path, '*.fastq')):
	os.system('perl prinseq-lite.pl -fastq ' + filename + ' -out_format 1 -out_good ' + filename.replace('.fastq',''))
	 