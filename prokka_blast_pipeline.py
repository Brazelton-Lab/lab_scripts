#! /usr/bin/env python

"""
BLAST sequences from PROKKA annotations

Takes a PROKKA FAA file and a file with an one or more IDs (each ID must be
on its own line). If a PROKKA annotation description contains an ID in
the IDs file, its sequence is BLASTed against NCBI BLAST and the results are
written to a custom TSV. This program is tailored to streamline part
of Katrina Twing's data analysis and the TSV is tailored to match.

Copyright:

    prokka_blast_pipeline  BLAST sequences from PROKKA annotations

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

from __future__ import division
import argparse
from bio_utils.iterators import fasta_iter, GFF3Reader
from Bio.Blast.NCBIWWW import qblast
from Bio.Blast import NCBIXML
from collections import defaultdict
import errno
import os
from socket import error as socket_error
import sys
from tqdm import tqdm

__author__ = 'Alex Hyer'
__email__ = 'theonehyer@gmail.com'
__license__ = 'GPLv3'
__maintainer__ = 'Alex Hyer'
__status__ = 'Production'
__version__ = '1.1.2'


def main(args):
    """Run program

    Args:
        args (NameSpace): ArgParse arguments dictating program use
    """

    tqdm.write('>>> Starting prokka_blast_pipeline')

    # Get Gene IDs from ID file
    ids = [gene_id.strip() for gene_id in args.id]

    tqdm.write('>>> Found {0} ID(s) in {1}'
               .format(str(len(ids)), args.id.name))

    # Get contigs that contain a feature ID matching a given Gene ID
    contigs = []
    gff_reader= GFF3Reader(args.gff3)
    for entry in gff_reader.iterate():
        try:  # Ignore features without a gene_feature field
            if entry.attributes['gene_feature'] in ids:
                contigs.append(entry.seqid)
        except KeyError:
            continue

    tqdm.write('>>> Found {0} contig(s) containing gene features matching '
               'given ID(s) in {1}'.format(str(len(contigs)), args.gff3.name))

    # Obtain PROKKA IDs and annotations of genes on contigs obtained earlier
    prokka_to_contig = defaultdict(str)
    prokka_to_gene = defaultdict(str)
    args.gff3.seek(0)
    for entry in gff_reader.iterate():
        if entry.seqid in contigs and 'gene_feature' in entry.attributes:
            prokka_to_contig[entry.attributes['ID']] = entry.seqid
            prokka_to_gene[entry.attributes['ID']] = entry.attributes[
                'gene_feature']

    tqdm.write('>>> Found {0} gene feature(s) on contigs matching given ID(s)'
               .format(str(len(prokka_to_contig))))

    # Get sequences from FAA file if they match a PROKKA ID obtained above
    blast_entries = []
    for entry in fasta_iter(args.faa):
        if entry.id in prokka_to_contig or ids[0] == '*':
            blast_entries.append(entry)

    tqdm.write('>>> Obtained {0} amino acid sequence(s) from {1}'
               .format(str(len(blast_entries)), args.faa.name))

    # Output header line
    args.output.write('Contig\tPROKKA_ID\tAnnotation\tGene ID\tSubject\t'
                      'Query Coverage(%)\tE-Value\tIdentity(%){0}'.format(
                                                                   os.linesep))
    tqdm.write('>>> BLASTing {0} amino acid sequence(s) against the NCBI {1} '
               'database'.format(str(len(blast_entries)), args.database))

    # BLAST sequences and calculate various summary values
    count = 0
    for entry in tqdm(blast_entries):
        # Continue attempting same blast until it succeeds w/o NCBI errors
        while True:
            try:
                result_handle = qblast(args.program, args.database,
                                       entry.sequence, alignments=args.top,
                                       descriptions=args.top,
                                       hitlist_size=args.top,
                                       expect=args.e_value)
                break
            except (ValueError, socket_error):  # Ignore NCBI
                continue

        # Process BLAST results
        result_generator = NCBIXML.parse(result_handle)
        for result in result_generator:
            for alignment in result.alignments:
                count += 1
                for hsp in alignment.hsps:
                    cov = float(hsp.align_length / len(entry.sequence)) * 100.0
                    perc = float(hsp.identities / len(entry.sequence)) * 100.0
                    taxonomy = alignment.hit_def.split('[')[1]
                    taxonomy = taxonomy.split(']')[0]

                    # Format and write output to custom TSV
                    output = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}{8}'\
                             .format(prokka_to_contig[entry.id], entry.id,
                                     entry.description,
                                     prokka_to_gene[entry.id], taxonomy,
                                     str(cov), str(hsp.expect), str(perc),
                                     os.linesep)
                    args.output.write(output)

    tqdm.write('>>> Wrote {0} total hit(s) to {1}'.format(str(count),
                                                          args.output.name))

    tqdm.write('>>> Exiting prokka_blast_pipeline')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.
                                     RawDescriptionHelpFormatter)
    parser.add_argument('-d', '--database',
                        default='nr',
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
                        help='FAA file from PROKKA containing amino acid '
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
                        default='blastp',
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
