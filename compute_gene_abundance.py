#!/usr/bin/env python
"""
Finds abundance of genes for metagenome on metagenome

Usage:

    compute_gene_abundance [--log_file] [--verbose] [--version] [--verify] 
                   [--fastr] [--fasta] [--normalize] <gff> <output>
"""

from __future__ import division

import argparse
from metameta.metameta_utils.output import output
import re
import statistics
import sys
from math import fabs

def gff3_to_dict(gff3_file, database):
    gff3_dict = {}
    with open(gff3_file, 'rU') as gff3_handle:
        for entry in gff3_iter(gff3_handle):
            db_id = extract_db_id(entry['attributes'], database)
            gff3_dict[entry['seqid']][db_id] = entry
    return gff3_dict

def fasta():
    output('Reading ' + args.gff[0], args.verbosity, 1,\
           log_file = args.log_file)
    gffFile = gff_dict(args.gff[0])
    toWrite = {}
    with open(args.output + '.tsv', 'w') as out_handle:
        for entry in entry_generator(args.fasta[0], 'fasta'):
            header = entry[0].lstrip('>').rstrip('\n').split()
            # catch cases when there is no read depth information
            entryId = header[0]
            readDepth = header[-1].split('_')[-1]
            if readDepth == '0':
                continue
            for hit in gffFile[entryId]:
                try:
                    gffDbIdSegment = re.findall('similar to AA sequence:.+;',\
                                                hit[-1])[0]
                    gffDbId = re.split('similar to AA sequence:|;', gffDbIdSegment)[1]
                    gffDbId = gffDbId.split(args.database + ':')[1]
                except IndexError:
                    continue
                else:
                    if args.normalize:
                        length = fabs(float(hit[4]) - float(hit[3]))
                        readDepth = normalize(readDepth, length)
                    toWrite[gffDbId] = toWrite.get(gffDbId, 0) + readDepth
                    
        for key in toWrite:
            out_handle.write(key + '\t' + str(toWrite[key]) + '\n')
                

def fastr():
    output('Reading ' + args.gff[0], args.verbosity, 1,\
           log_file = args.log_file)
    gffFile = gff_dict(args.gff[0])
    with open(args.output + '.tsv', 'w') as out_handle:
        for entry in entry_generator(args.fastr[0], 'fastr'):
            header = entry[0].lstrip('+').rstrip('\n')
            if header in gffFile:
                for ann in gffFile[header]:
                    gffStart = int(ann[3]) - 1
                    gffEnd = int(ann[4]) - 1
                    depthData = decompress_fastr(entry[1]).split('-')
                    depthData = [int(i) for i in depthData[gffStart:gffEnd]]
                    if len(depthData) > 0:
                        average = statistics.mean(depthData)
                        if args.normalize:
                            length = len(depthData)
                            average = normalize(average, length)
                        try:
                            gffDbIdSegment = re.findall('similar to AA sequence:.+;',\
                                                        ann[-1])[0]
                            gffDbId = re.split('similar to AA sequence:|;', gffDbIdSegment)[1]
                            gffDbId = gffDbId.split(args.database + ':')[1]
                            try:
                                assert float(gffDbId)
                            except:
                                if gffDbId in toWrite:
                                    if toWrite[gffDbId] < average:
                                        toWrite[gffDbId] = average
                                else:
                                    toWrite[gffDbId] = average
                        except IndexError:
                            pass
        for key in toWrite:
            out_handle.write(key + '\t' + str(toWrite[key]) + '\n')

def normalize(read_depth, length):
    kb = float(length)/1000.0
    rpk = float(read_depth)/kb
    return rpk

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = __doc__,
                                        formatter_class = argparse.\
                                        RawDescriptionHelpFormatter)
    parser.add_argument('gff',
                        default = None,
                        nargs = '?',
                        help = 'GFF3 file with annotations for the same'\
                        + ' metagenome as the FASTR file')
    parser.add_argument('output', metavar='out_prefix',
                        default = None,
                        nargs = '?',
                        help = 'name of GFF3 file to write')
    parser.add_argument('--fasta', metavar='FASTA',
                        default = None,
                        nargs = '?',
                        help = 'IDBA generated FASTA file')
    parser.add_argument('--fastr', metavar='FASTR',
                        default = None,
                        nargs = '?',
                        help = 'FASTR file containing read depth data'\
                        + ' for a metatranscriptome mapped onto a metagenome')
    parser.add_argument('--normalize',
                        default = None,
                        action = 'store_true',
                        help = 'normalize the data with afk')
    parser.add_argument('-l', '--log_file', metavar='LOG',
                        default = None,
                        help = 'log file to print all messages to')
    parser.add_argument('-v', '--verbosity',
                        action = 'count',
                        default = 0,
                        help = 'increase output verbosity')
    parser.add_argument('--verify',
                        action = 'store_true',
                        help = 'verify input files before use')
    parser.add_argument('--version',
                        action = 'store_true',
                        help = 'prints tool version and exits')
    parser.add_argument('-d', '--database', metavar='DB',
                        type=str,
                        default='ko',
                        help="extract hits from database \'DB\' [default: ko]")
    args = parser.parse_args()

    if args.version:
        print(__version__)
    elif args.fastr == None and args.fasta == None:
        print(__doc__)
    elif args.gff == None:
        message = 'Must specify a FASTA/R, GFF3, and output file'
        output(message, args.verbosity, 0,\
               log_file = args.log_file)
    else:
        if args.fastr != None:
            fastr()
        if args.fasta != None:
            fasta()

    sys.exit(0)
