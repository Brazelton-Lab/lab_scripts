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
