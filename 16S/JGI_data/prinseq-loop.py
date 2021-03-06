#! /usr/bin/env python

"""

Copyright:

    prinseq-loop  runs prinseq-lite on all .fastq files in a folder

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

path = r'./'
for filename in glob.glob(os.path.join(path, '*.fastq')):
	os.system('prinseq-lite -fastq ' + filename + ' -out_format 1 -out_good ' + filename.replace('.fastq',''))
	 
