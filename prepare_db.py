#! /usr/bin/env python
"""
For internal preparation of reference databases. Output depends on the given
reference source.

Copyright:

    prepare_db  For internal preparation of reference databases

    Copyright (C) 2016  William Brazelton, Christopher Thornton

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

import argparse
from bio_utils.iterators import fasta_iter
from bz2 import BZ2File
from gzip import GzipFile
from lzma import LZMAFile
import json
import io
import os
import re
import sys

__author__ = "Christopher Thornton"
__license__ = 'GPLv2'
__maintainer__ = 'Christopher Thornton'
__status__ = "Beta"
__version__ = "0.1.0"


class Open(argparse.Action):
    """Argparse Action that detects and opens compressed files for rw

    Attributes:
        option_strings (list): list of str giving command line flags that
                               call this action

        dest (str): Namespace reference to value

        mode (str): mode to pass to (de)compression algorithm

        nargs (bool): True if multiple arguments specified

        **kwargs (various): optional arguments to pass to argparse and algo
    """

    def __init__(self, option_strings, dest, mode='rb', nargs=None, **kwargs):
        """Initialize class and spawn self as Base Class w/o nargs

        Warns:
            ImportError: if Open cannot import a compression library,
                         it warns the user that it cannot open the
                         corresponding file type

        Raises:
            ValueError: if nargs is not None, Open does not accept nargs
        """
        from importlib import import_module
        from warnings import warn
        from os import linesep

        # Only accept a single value to analyze
        if nargs is not None:
           raise ValueError('nargs not allowed for Open')

        # Call self again but without nargs
        super(Open, self).__init__(option_strings, dest, **kwargs)

        # Store and establish variables used in __call__
        self.kwargs = kwargs
        self.mode = mode.lower().strip()
        self.modules = {}

        modules_to_import = {
            'bz2': 'BZ2File',
            'gzip': 'GzipFile',
            'lzma': 'LZMAFile'
        }

        # Dynamically import compression libraries and warn about failures
        for mod, _class in modules_to_import.items():
            try:
                self.modules[_class] = getattr(import_module(mod), _class)
            except (ImportError, AttributeError) as e:
                self.modules[_class] = open
                warn('Cannot process {0} files due to following error:'
                     '{1}{2}{1}You will need to install the {0} library to '
                     'properly use these files. Currently, such files will '
                     'open in text mode.'.format(mod, linesep, e))
    # Credits: https://stackoverflow.com/questions/13044562/
    # python-mechanism-to-identify-compressed-file-type-and-uncompress
    def __call__(self, parser, namespace, value, option_string=None, **kwargs):
        """Detects and opens compressed files

        Args:
            parser (ArgumentParser): parser used to generate values

            namespace (Namespace): namespace to set values for

            value (str): actual value specified by user

            option_string (str): argument flag used to call this function

            **kwargs (various): optional arguments later passed to the
                                compression algorithm
        """
        from inspect import getfullargspec
        import io

        filename = value  # For readability

        algo = io.open  # Default to plaintext

        algo_map = {
            'bz2': self.modules['BZ2File'],
            'gz':  self.modules['GzipFile'],
            'xz':  self.modules['LZMAFile']
        }

        # Base compression algorithm on file extension
        ext = value.split('.')[-1]
        try:
            algo = algo_map[ext]
        except KeyError:
            pass

        # Filter all **kwargs by the args accepted by the compression algo
        algo_args = set(getfullargspec(algo).args)
        good_args = set(self.kwargs.keys()).intersection(algo_args)
        _kwargs = {arg: self.kwargs[arg] for arg in good_args}


        # Open the file using parameters defined above and store in namespace
        try:
            handle = algo(value, mode=self.mode, **_kwargs)
        except ValueError:
            mode = self.mode.lstrip('U')[0]
            handle = io.TextIOWrapper(algo(value, mode=mode, **_kwargs), encoding='utf-8')

        setattr(namespace, self.dest, handle)


def open_io(infile, mode='rb'):
    """
    Open input or output files, using the file extension to dictate how the 
    file should be opened
    """
    algo = io.open  # Default to plaintext

    algo_map = {
        'bz2': BZ2File,
        'gz':  GzipFile,
        'xz':  LZMAFile
    }

    # Base compression algorithm on file extension
    ext = infile.split('.')[-1]
    try:
        algo = algo_map[ext]
    except KeyError:
        pass

    handle = io.TextIOWrapper(algo(infile, mode=mode), encoding='utf-8')

    return handle


def sub_foam(args, line=None):
    """
    Subcommand for generating internal reference files for the FOAM database
    """
    in_h = args.foam_in
    out_h = args.out_map

    r = re.compile("((?<=\d\d_).+)")

    append = list.append
    join = str.join
    strip = str.strip

    if line is None:
        line = next(in_h)  # Read header

    # Check if input is text or bytestream
    if (isinstance(line, bytes)):
        def next_line(i):
            return next(i).decode('utf-8')
    else:
        next_line = next

    current_path = ''
    accessions = []
    try:
        while True:  #loop until StopIteration raised

            line = strip(next_line(in_h)).split('\t')

            # pathway levels are first four columns with corresponding KO in fifth
            path = line[:-1]
            acc = line[-1]

            # remove L1 FOAM entry numbers from pathway
            try:
                # find to the lowest resolution available for a given pathway
                model = r.search(join(';', path[: path.index('')])).group()
            except ValueError:
                # pathway has values for all four FOAM levels
                model = r.search(join(';', path)).group()

            if model == current_path:
                append(accessions, acc)
                continue
        
            # output pathway map of previous path
            if accessions:
                output = "{}\t{}\n".format(current_path, '\t'.join(accessions))
                out_h.write(output)

            current_path = model  #new path is now current path

            # start list of accessions for current path
            accessions = [acc]

    except StopIteration:  #output last pathway
        output = "{}\t{}\n".format(model, '\t'.join(accessions))
        out_h.write(output)


def sub_card(args):
    """
    Subcommand for generating internal reference files for the CARD database
    """

    in_h = args.card_in
    out_h = args.out_map

    supported_models = args.card_models

    db_version = ' v{}'.format(args.db_version) if args.db_version else ''
    ref_db = "CARD{}".format(db_version)

    protein_models = ["protein homolog model", "protein variant model", \
                      "protein overexpression model", "protein knockout model"]
    meta_data = {}

    out_map = {}
    for mtype in supported_models:
        out_fna = open("CARD-{}.fna".format(mtype.replace(' ', '_')), 'wt') \
                  if args.fasta else ""
        out_faa = open("CARD-{}.faa".format(mtype.replace(' ', '_')), 'wt') \
                  if args.fasta and mtype in protein_models else ""
        out_map[mtype] = {"fna": out_fna, "faa": out_faa, "counts": 0}

    # Speed tricks
    append = list.append
    join = str.join

    # Load CARD database
    json_db = json.load(in_h)

    db_totals = 0
    ontologies = {}
    for model in json_db:
        try:
            aro = json_db[model]["ARO_accession"]
        except KeyError:
            continue
        except TypeError:
            continue

        db_totals += 1

        model_id = json_db[model]["model_id"]
        model_name = json_db[model]["model_name"]
        product = json_db[model]["ARO_description"]
        mechanisms = []
        drugs = []
        drug_classes = []
        snps = []

        model_type = json_db[model]["model_type"]
        if model_type not in supported_models:
            continue
        else:
            out_map[model_type]["counts"] += 1

        # Parameter values
        try:
            params = json_db[model]["model_param"]
        except KeyError:
            continue

        if "blastp_bit_score" in params and "blastn_bit_score" in params:
            print("model ID {} has both blastp and blastn bitscore values".format(model_id))
            sys.exit(1)

        try:
            scores = params["blastp_bit_score"]["param_value"]
        except KeyError:
            try:
                scores = params["blastn_bit_score"]["param_value"]
            except KeyError:
                scores = ""

        if "snp" in params:
            snpkeys = params["snp"]["param_value"]
            for snpkey in snpkeys:
                snps.append(snpkeys[snpkey])

        # CARD Categories
        categories = json_db[model]["ARO_category"]
        for category in categories:
            category_name = categories[category]["category_aro_name"]

            class_name = categories[category]["category_aro_class_name"]

            if class_name == "Resistance Mechanism":
                mechanisms.append(category_name)
            elif class_name == "Drug Class":
                drug_classes.append(category_name)
            elif class_name == "Antibiotic":
                drugs.append(category_name)

        # CARD sequences
        seqs = json_db[model]["model_sequences"]["sequence"]
        for record in seqs:
            acc = "{}_{}".format(model_id, record)
            organism = seqs[record]["NCBI_taxonomy"]["NCBI_taxonomy_name"]
            nucl_seq = seqs[record]["dna_sequence"]["sequence"]
            seq_acc_nucl = seqs[record]["dna_sequence"]["accession"]
            seqlen = len(nucl_seq)

            if model_type in protein_models:
                prot_seq = seqs[record]["protein_sequence"]["sequence"]

            if acc not in meta_data:
                meta_data[acc] = {
                                  'organism': organism,
                                  'product': product,
                                  'gene_length': seqlen,
                                  'gene': model_name,
                                  'accession': seq_acc_nucl,
                                  'gene_family': aro,
                                  'database': ref_db,
                                  'bitscore': scores,
                                  'snp': snps,
                                  'resistance_mechanism': mechanisms,
                                  'compound': drugs,
                                  'drug_class': drug_classes
                                  }
            else:
                print("error: redundant model {}".format(acc), file=sys.stderr)
                sys.exit(1)

            # Generate CARD fasta files by model type, if requested
            if args.fasta:

                try:
                    model_map = out_map[model_type]
                except KeyError:
                    continue

                out_fna = model_map["fna"]
                out_faa = model_map["faa"]

                if out_fna:
                    out_fna.write(">{}\n{}\n".format(acc, nucl_seq))

                if out_faa:
                    out_faa.write(">{}\n{}\n".format(acc, prot_seq))

    # Output internal relational database
    json.dump(meta_data, out_h, sort_keys=True, indent=4, separators=(',', ': '))

    # Database statistics:
    print("Total database size: {}".format(db_totals), file=sys.stderr)
    for model_type in out_map:
        print("  - models in '{}':\t{!s}".format(model_type, out_map[model_type]["counts"]), file=sys.stderr)


def sub_kegg(args):
    """
    Subcommand for generating internal reference files for the KEGG database
    """

    out_h = args.out_map
    taxonomy = parse_kegg_taxonomy(args.kegg_tax) if args.kegg_tax else None

    db_version = ' v{}'.format(args.db_version) if args.db_version else ''
    ref_db = "KEGG{}".format(db_version)

    kegg_map = {}  # store KEGG entries

    # Parse DAT files, if provided
    if args.kegg_dat:

        dat_totals = 0
        for line in dat_h:
            split_line = line.split('\t')
            try:
                acc, ko, aa_len, product = split_line
            except ValueError:
                num_col = len(split_line)
                print("error: unknown file format for {}. Expected "
                      "4 columns, but only {} were provided"
                      .format(dat_file, num_col), file=sys.stderr)
                sys.exit(1)

            dat_totals += 1

            tax_code = acc.split(':')[0]
            try:
                organism = taxonomy[tax_code]
            except TypeError:  #no taxonomy file provided
                organism = ''
            except KeyError:  #no entry in taxonomy file for tax_code
                print("error: code {} from {} not found in taxonomy file"\
                      .format(tax_code, acc), file=sys.stderr)
                organism = ''

            kegg_map[acc] = {
                             'organism': organism,
                             'product': '',
                             'gene_length': int(aa_len) * 3,
                             'gene': '',
                             'gene_family': ko,
                             'database': ref_db,
                             }

    # Parse FASTA files
    fasta_totals = 0
    for record in fasta_iter(fasta_h):
        fasta_totals += 1

        acc = record.id

        split_desc = record.description.split(';')
        try:
            gene, product = split_desc
        except ValueError:
            gene = ''
            product = record.description
        else:
            product = product.rstrip()

        gene_len = len(record.sequence) * 3

        if acc in kegg_map:
            if kegg_map[acc]['gene_length'] != gene_len:
                print("error: accession {} has differing gene lengths in "
                      "FASTA and DAT files".format(acc), file=sys.stderr)
                sys.exit(1)

            kegg_map[acc]['product'] = product
            kegg_map[acc]['gene'] = gene
        else:
            if args.kegg_dat:
                print("error: accession {} in FASTA file does not have a "
                      "corresponding entry in DAT file".format(acc), \
                      file=sys.stderr)
                sys.exit(1)
            else:
                kegg_map[acc] = {
                                 'organism': taxonomy[acc],
                                 'product': product,
                                 'gene_length': gene_len,
                                 'gene': gene,
                                 'gene_family': '',
                                 'database': ref_db,
                                 }

    if args.kegg_dat:
        if dat_totals != fasta_totals:
            print("error: DAT file and FASTA file should contain the same number "
                  "of entries", file=sys.stderr)
            sys.exit(1)

    # Output internal relational database
    json.dump(kegg_map, out_h, sort_keys=True, indent=4, separators=(',', ': '))


def sub_uniprot(args):
    """
    Subcommand for generating internal reference files for the UniProt database
    """
    pass


def sub_bacmet(args):
    """
    Subcommand for generating internal reference files for the BacMet database
    """
    PRODUCT = re.compile('.+?(?=OS\=)')
    ORG = re.compile('(?<=OS\=).+?(?=[A-Z]{2}\=)')

    in_h = args.bacmet_in
    in_m = args.bacmet_map
    out_h = args.out_map
    out_f = args.bacmet_fa

    db_version = ' v{}'.format(args.db_version) if args.db_version else ''
    ref_db = "BacMet{}".format(db_version)

    # Parse BacMet mapping info
    meta_data = {}
    for line in in_m:
        line = line.decode('utf-8')

        if line.startswith('#'):
            continue

        split_line = line.strip().split('\t')
        try:
            acc, ident, gene, location, organism, drugs, desc = split_line
        except ValueError:
            num_cols = len(split_line)
            print("error: unknown mapping file format. Seven columns expected, "
                  "{!s} columns provided".format(num_cols), file=sys.stderr)
            sys.exit(1)

        drugs = drugs.split(',')
        drugs = [i.lstrip() for i in drugs]

        meta_data[ident] = {'accession': acc,
                           'gene': gene,
                           'gene_length': '',
                           'database': ref_db,
                           'organism': organism,
                           'product': '',
                           'compound': drugs
                          }

    # Parse BacMet FASTA file
    for record in fasta_iter(in_h):
        header = record.id.split('|')
        ident = header[0]
        seq_len = len(record.sequence)

        try:
            product = PRODUCT.search(record.description).group()
        except AttributeError:
            product = record.description
        product = product.strip()

        if ident in meta_data:
            meta_data[ident]['gene_length'] = seq_len * 3
            meta_data[ident]['product'] = product
        else:
            try:
                organism = ORG.search(record.description).group()
            except AttributeError:
                organism = ''
            organism = organism.strip()

            gene = header[1]
            acc = header[3]
            meta_data[ident] = {'accession': acc,
                                'gene': gene,
                                'gene_length': seq_len * 3,
                                'database': ref_db,
                                'organism': organism,
                                'product': product,
                                'compound': []
                               }

        # Output FASTA with simplified headers
        out_f.write(">{}\n{}\n".format(ident, record.sequence))

    # Output internal relational database
    json.dump(meta_data, out_h, sort_keys=True, indent=4, separators=(',', ': '))


def reaction_alts(acc_list):
    """
    Convert list of proteins using standard convention for alternative 
    enzymes/reactions
    """
    if len(acc_list) > 1:
        return("({})".format(",".join(acc_list)))
    else:
        return(acc_list[0])


def do_nothing(*args):
    pass


def parse_kegg_taxonomy(in_h):
    tax = {}
    for line in in_h:
        line = line.decode('utf-8')

        if line.startswith('#'):
            continue

        split_line = line.strip().split('\t')
        try:
            ident_1, tax_code, ident_2, taxonomy = split_line
        except ValueError:
            num_cols = len(split_line)
            print("error: unknown taxonomy file format. Four columns expected, "
                  "{!s} columns provided".format(num_cols), file=sys.stderr)
            sys.exit(1)

        tax[tax_code] = taxonomy

    return tax


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parent_parser = argparse.ArgumentParser(add_help=False)
    subparsers = parser.add_subparsers(title="subcommands",
        help="exactly one of these commands is "
             "required")
    parent_parser.add_argument('-v', '--db-version',
        type=str,
        help="database version")
    parent_parser.add_argument('-o', '--out',
        metavar='OUTFILE',
        dest='out_map',
        action=Open,
        mode='wt',
        default=sys.stdout,
        help="output relational database [default: output to stdout]")
    # FOAM-specific arguments
    foam_parser = subparsers.add_parser('foam',
        parents=[parent_parser],
        help="generate internal files for the FOAM ontology")
    foam_args = foam_parser.add_argument_group("FOAM-specific arguments")
    foam_args.add_argument('-f', '--foam',
        metavar='INFILE',
        dest='foam_in',
        action=Open,
        mode='rb',
        help="input tabular FOAM ontology file relating Kegg Orthologs (KOs) "
             "to a hierarchical representation of biochemical functions and "
             "pathways")
    foam_parser.set_defaults(func=sub_foam)
    # CARD-specific arguments
    card_parser = subparsers.add_parser('card',
        parents=[parent_parser],
        help="generate internal files for the CARD reference database. CARD "
             "is not as established a database as KEGG or UniProt. Future "
             "updates may break compatibility, requiring modifications to the "
             "program")
    card_args = card_parser.add_argument_group("CARD-specific arguments")
    card_args.add_argument('-c', '--card',
        metavar='card.json',
        dest='card_in',
        action=Open,
        mode='rb',
        help="input complete CARD database in JSON format")
    card_args.add_argument('-cm', '--card-models',
        metavar='MODEL [,MODEL]',
        dest='card_models',
        default=["protein homolog model", "protein variant model", \
                 "rRNA gene variant model", "protein overexpression model", \
                 "protein knockout model"],
        help="Output information on the following CARD model types [default: "
              "'protein homolog model', 'protein variant model', 'rRNA gene "
              "variant model', 'protein overexpression model', 'protein "
              "knockout model']")
    card_args.add_argument('--card-fasta',
        dest='fasta',
        action='store_true',
        help="generate FASTA files of reference sequences by model type")
    card_parser.set_defaults(func=sub_card)
    # KEGG-specific arguments
    kegg_parser = subparsers.add_parser('kegg',
        parents=[parent_parser],
        help="generate internal files for KEGG")
    kegg_args = kegg_parser.add_argument_group("KEGG-specific arguments")
    kegg_args.add_argument('-k', '--kegg',
        metavar='INFILE',
        dest='kegg_fasta',
        action=Open,
        mode='rb',
        help="input FASTA file of KEGG protein sequences")
    kegg_args.add_argument('-kd', '--kegg-dat',
        metavar='INFILE',
        dest='kegg_dat',
        action=Open,
        mode='rb',
        help="input DBGET flat file of KEGG genes mapped to their KO "
             "assignments")
    kegg_args.add_argument('-kt', '--kegg-taxonomy',
        metavar='INFILE',
        dest='kegg_tax',
        action=Open,
        mode='rb',
        help="input tabular taxonomy file")
    kegg_parser.set_defaults(func=sub_kegg)
    # UniProt-specific arguments
    uniprot_parser = subparsers.add_parser('uniprot',
        parents=[parent_parser],
        help="generate internal files for the UniProt protein databases")
    uniprot_args = uniprot_parser.add_argument_group("UniProt-specific arguments")
    uniprot_args.add_argument('-u', '--uniprot',
        metavar='INFILE',
        dest='uniprot_in',
        action=Open,
        mode='rb',
        help="input DBGET flat file containing UniProt entries")
    uniprot_parser.set_defaults(func=sub_uniprot)
    # BacMet-specific arguments
    bacmet_parser = subparsers.add_parser('bacmet',
        parents=[parent_parser],
        help="generate internal files for the BacMet reference databases")
    bacmet_args = bacmet_parser.add_argument_group("BacMet-specific arguments")
    bacmet_args.add_argument('-b', '--bacmet',
        metavar='INFILE',
        dest='bacmet_in',
        action=Open,
        mode='rb',
        help="input FASTA file of BacMet protein sequences")
    bacmet_args.add_argument('-bm', '--bacmet-map',
        metavar='MAP',
        dest='bacmet_map',
        action=Open,
        mode='rb',
        help="input BacMet mapping file")
    bacmet_args.add_argument('-bf', '--bacmet-fasta',
        dest='bacmet_fa',
        action=Open,
        mode='wt',
        default="BacMet.fa",
        help="generate FASTA file with simplified headers [default: "
             "./BacMet.fa]")
    bacmet_parser.set_defaults(func=sub_bacmet)
    args = parser.parse_args()

    args.func(args)


if __name__ == "__main__":
    main()
    sys.exit(0)
