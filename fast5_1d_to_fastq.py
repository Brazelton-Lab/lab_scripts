#! /usr/bin/env python

"""
Copyright:

    fast5_1d_to_fasta.py Extract 1D read data from PORESEQ FAST5 file
    Copyright (C) 2016  William Brazelton, Alex Hyer

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

from __future__ import print_function
import h5py
from Bio import SeqIO
from StringIO import StringIO
import sys

for file in sys.argv[1:]:
    hdf = h5py.File(file, 'r')
    try:
        fq = hdf['Analyses']['Basecall_1D_CDNA_000']['BaseCalled_template']['Fastq'][()]
        print(fq.strip())
    except Exception, e:
        pass
    hdf.close()
