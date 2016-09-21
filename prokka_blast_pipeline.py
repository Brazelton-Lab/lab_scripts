#! /usr/bin/env python

from __future__ import division
from __future__ import print_function

"""BLAST sequences from PROKKA annotations

Takes a PROKKA FNA file and a file with an one or more IDs (each ID must be
on its own line). If a PROKKA annotation description contains an ID in
the IDs file, its sequence is BLASTed against NCBI BLAST and the results are
written to a custom TSV.
"""

import argparse
from bio_utils.iterators import fasta_iter, gff3_iter
from Bio.Blast.NCBIWWW import qblast
from Bio.Blast import NCBIXML
import os
import sys

__author__ = 'Alex Hyer'
__email__ = 'theonehyer@gmail.com'
__license__ = 'GPLv3'
__maintainer__ = 'Alex Hyer'
__status__ = 'Alpha'
__version__ = '0.0.1a4'


def main(args):
    """Run program

    Args:
        args (NameSpace): ArgParse arguments dictation program use
    """

    # Get IDs from ID file
    ids = [gene_id for gene_id in args.id]

    # Get PROKKA Ids a feature ID matches a given Gene ID
    prokka_ids = []
    for entry in gff3_iter(args.gff3):
        if entry.attributes['gene_feature'] in ids and \
                entry.attributes['ID'] not in prokka_ids:
            prokka_ids.append(entry.attributes['ID'])

    # Get sequences from FAA file if they match a PROKKA ID
    entries = []
    for entry in fasta_iter(args.faa):
        if entry.id in prokka_ids or gene_id == '*':
            if entry.id not in [seq[0].id for seq in entries]:
                entries.append((entry, gene_id))

    # Output header line
    args.output.write('Contig\tPROKKA_ID\tAnnotation\tGene ID\tSubject\t'
                      'Query Coverage\tE-Value\tIdentity{0}'.format(
                                                                   os.linesep))

    # BLAST sequences
    for entry in entries:
        result_handle = qblast(args.program, args.database, entry[0].sequence,
                               alignments=args.top, descriptions=args.top,
                               hitlist_size=args.top, expect=args.e_value)
        for alignment in NCBIXML.parse(result_handle).alignments:
            for hsp in alignment.hsps:
                prokka_id = entry.description.split(' ')[1]
                ann = entry.description.split(' ')[2]
                cov = float(hsp.align_length / len(entry[0].sequence)) / 100.0
                output = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}{8}'.format(
                    entry.id, prokka_id, ann, entry[1], hsp.sbjct, str(cov),
                    str(hsp.expect), str(hsp.identities), os.linesep)
                args.output.write(output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.
                                     RawDescriptionHelpFormatter)
    parser.add_argument('-d', '--database',
                        default='nt',
                        choices=[
                            'nr',
                            'refseq_protein',
                            'swissprot',
                            'env_nr',
                            'tsa_nr'
                        ],
                        help='BLAST database to use. See Table 2 from '
                             'ftp://ftp.ncbi.nlm.nih.gov/pub/factsheets/HowTo_BLASTGuide.pdf')
    parser.add_argument('-e', '--e_value',
                        default=10.0,
                        type=float,
                        help='Maximum E-Value of alignment permitted')
    parser.add_argument('-f', '--faa',
                        required=True,
                        type=argparse.FileType('r'),
                        help='FAA file from PROKKA containing nucleotide '
                             'sequences of annotated proteins')
    parser.add_argument('-g', '--gff3',
                        required=True,
                        type=argparse.FileType('r'),
                        help='GFF3 file from PROKKA accompanying the FAA file')
    parser.add_argument('-i', '--id',
                        required=True,
                        type=argparse.FileType('r'),
                        help='line-separated list of IDs designating which '
                             'PROKKA IDs to BLAST with. If only ID is "*" '
                             'will match all IDs.')
    parser.add_argument('-o', '--output',
                        type=argparse.FileType('w'),
                        help='file name for output tsv')
    parser.add_argument('-p', '--program',
                        default='blastn',
                        type=str,
                        choices=[
                            'blastp',
                            'tblastn'
                        ],
                        help='BLAST+ program to search with')
    parser.add_argument('-t', '--top',
                        default=1,
                        type=int,
                        help='number of best BLAST results to return')
    args = parser.parse_args()

    main(args)

    sys.exit(0)
