#! /usr/bin/env python3

"""Convert various file types to TSVs for use with anvi'o

This program has essentially two functions: splitting GFF3 files
into a gene_location and genes file to tell anvi'o both where your gene calls
are and what they consist of, and to take FASTA files and produce a tsv for
anvi'o describing which contig is in which bin.

Updated in 2020 to handle prodigal GFF3 files and their annotations, instead of prokka GFF3 files that the previous version of this program was customized for.

Copyright:

	anvi_converter.py Convert various file types to TSVs for use with anvi'o
	Copyright (C) 2020	William Brazelton, Alex Hyer

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

import argparse
from bio_utils.iterators import fasta_iter, GFF3Reader
import os
import sys

__author__ = 'Alex Hyer'
__email__ = 'theonehyer@gmail.com'
__license__ = 'GPLv3'
__maintainer__ = 'Alex Hyer'
__status__ = 'Production'
__version__ = '1.1.10'


def main(args):
	"""Run program

	Args:
		 args (NameSpace): ArgParse arguments controlling program flow
	"""

	if args.tool == 'bins':
		with open(args.bins) as fasta: 				# changed this code to accept a file containing a list of names
			for line in fasta:
				filename = args.path2bins + '/' + line.strip()
				with open(filename, 'r') as file_handle:
					for entry in fasta_iter(file_handle):
						file_name = os.path.basename(file_handle.name)
						file_name = '.'.join(file_name.split('.')[:-1])
						args.output.write('{0}\t{1}{2}'.format(entry.id,
															   file_name,
															   os.linesep))

	if args.tool == 'gff':

		gff_reader = GFF3Reader(args.gff3)

		naughty_gene_calls = {}
		for entry in fasta_iter(args.faa):
			if '*' in entry.sequence:
				naughty_gene_calls[entry.id] = ''

		with open(args.prefix + '.gene_locations.tsv', 'w') as lh, \
				open(args.prefix + '.genes.tsv', 'w') as gh:

			caller_id = 1
			lh.write('gene_callers_id\tcontig\tstart\tstop\t'
					 'direction\tpartial\tsource\tversion{0}'
					 .format(os.linesep))
			gh.write('gene_callers_id\tsource\taccession\t'
					 'function\te_value{0}'.format(os.linesep))

			for entry in gff_reader.iterate():

				if 'ID' not in entry.attributes.keys():
					continue

				if entry.attributes['ID'] in naughty_gene_calls:
					continue

				if entry.type == 'CDS': #and \				commented this out because I don't want to skip entries without a "gene" attribute
						#'gene' in entry.attributes.keys():

					# Reformat data for gene locations file
					direction = 'f' if entry.strand == '+' else 'r'
					program, version = entry.source.split('_')				# changed this from ":"
					lh.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\tv{7}{8}'
							 .format(str(caller_id),
									 entry.seqid,
									 str(entry.start - 1),
									 str(entry.end - 3),
									 direction,
									 '0', program, version, os.linesep))

					# Reformat data for genes file
					if 'Alias' in entry.attributes.keys():		# only proceed if specified field is present in attributes (but still want such genes to appear in gene_locations file above)
						gh.write('{0}\t{1}\t{2}\t{3}\t{4}{5}'
								 .format(str(caller_id),
										 entry.source,
										 entry.attributes['Alias'],			# changed this from 'gene' because the Alias field contains the actual KEGG seq name
										 entry.attributes['product'],
										 '0', os.linesep))

					caller_id += 1


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description=__doc__,
									 formatter_class=argparse.
									 RawDescriptionHelpFormatter)
	subparsers = parser.add_subparsers(title='Tool',
									   dest='tool')

	bins = subparsers.add_parser('bins',
								 help='read multiple FASTA files and produce '
									  'TSV relating entries to file name')
	bins.add_argument('-b','--bins',
					  type=str,
					  help='file containing a list of FASTA files where each file contains the contigs for one bin')
	bins.add_argument('-p','--path2bins',
					  type=str,
					  help='path to bins listed in provided file')
	bins.add_argument('-o','--output',
					  type=argparse.FileType('w'),
					  help='output file')

	gff = subparsers.add_parser('gff',
								   help='convert GFF3 file into two TSVs containing called gene locations and their annotations')
	gff.add_argument('-g','--gff3',
						type=argparse.FileType('r'),
						help='GFF3 file to convert')
	gff.add_argument('-f','--faa',
						type=argparse.FileType('r'),
						help='FAA file corresponding to GFF3 file')
	gff.add_argument('-x','--prefix',
						type=str,
						help='prefix for output files')
	args = parser.parse_args()

	main(args)

	sys.exit(0)
