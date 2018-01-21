#! /usr/bin/env python
"""
For internal preparation of reference databases. Output depends on the given
reference source.
"""

from __future__ import print_function

import argparse
from bio_utils.iterators import fasta_iter
from bz2 import BZ2File
from gzip import GzipFile
import json
import io
import os
import re
import sys

__author__ = "Christopher Thornton"
__license__ = 'GPLv2'
__maintainer__ = 'Christopher Thornton'
__status__ = "Alpha"
__version__ = "0.0.1"


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


def sub_foam(args, line=None):
    """
    Subcommand for generating internal reference files for the FOAM database
    """
    in_h = args.foam_in
    out_h = args.foam_out

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
    snp_re = re.compile("([A-Z])([0-9]+)([A-Z])")

    in_h = args.card_in
    out_h = args.card_out
    supported_models = args.card_models
    db_version = "1.2.1"

    if args.card_genes:
        meta_data = {}

    if args.card_fasta:
        out_map = {}
        for model_type in supported_models:
            outfile = "CARD-{}.fa".format(model_type.replace(' ', '_'))
            out_map[model_type] = {"outfile": open(outfile, 'wt'), "counts": 0}

    if args.card_snp:
        args.card_snp.write("#Accession\tPosition\tWild Type\tMutant\n")

    # Speed tricks
    append = list.append
    join = str.join

    # Load CARD database
    json_db = json.load(in_h)

    db_totals = 0
    ontologies = {}
    for model in json_db:
        ontology = []

        try:
            categories = json_db[model]["ARO_category"]
        except KeyError:  #handle database description
            continue
        except TypeError:  #handle version information
            continue

        model_type = json_db[model]["model_type"]
        if model_type not in supported_models:
            continue
        else:
            out_map[model_type]["counts"] += 1

        for category in categories:
            cat_name = categories[category]["category_aro_name"]
            append(ontology, cat_name)

        try:
            aro_name = json_db[model]["ARO_name"]
        except KeyError:
            print("No ARO name for:\n{}".format(json.dumps(json_db[model], \
                  sort_keys=True, indent=4, separators=(',', ': ')), \
                  file=sys.stderr))
            sys.exit(1)

        # No clear hierarchy supplied by CARD for ARO categories, so define 
        # one and sort categories by new hierarchy
        ontology = sort_by_level(ontology)
        ontology = join(';', ontology)

        try:
            acc = "ARO:{}".format(json_db[model]["ARO_accession"])
        except KeyError:
            print("No Accession for:\n{}".format(json.dumps(json_db[model], \
                  sort_keys=True, indent=4, separators=(',', ': ')), \
                  file=sys.stderr))
            sys.exit(1)
        else:
            db_totals += 1
        
        try:
            append(ontologies[ontology], acc)
        except KeyError:
            ontologies[ontology] = [acc]

        try:
            params = json_db[model]["model_param"]
        except KeyError:
            print("No model parameters available for {}".format(json.dumps(json_db[model], \
                  sort_keys=True, indent=4, separators=(',', ': '))), file=sys.stderr)
            sys.exit(1)

        try:
            sequences = json_db[model]["model_sequences"]
        except KeyError:
            print("No model sequence for {}".format(json.dumps(json_db[model], \
                  sort_keys=True, indent=4, separators=(',', ': '))), file=sys.stderr)
            sys.exit(1)

        entries = sequences["sequence"]
        for entry in entries:
            gene_id = entries[entry]["protein_sequence"]["accession"] if \
                      "protein" in model_type else \
                      entries[entry]["dna_sequence"]["accession"]

            # Generate SNPs file, if requested
            if args.card_snp:
                if "snp" in params:
                    snpkeys = params["snp"]["param_value"]
                    for snpkey in snpkeys:
                        wildtype, snp_pos, mutant = snp_re.search(snpkeys[snpkey]).groups()
                        args.card_snp.write("{}\t{}\t{}\t{}\n".format(gene_id, snp_pos, wildtype, mutant))

            # Generate CARD gene metadata
            if args.card_genes:
                # Alignment score thresholds
                if "blastp_bit_score" in params:
                    bitscore = params["blastp_bit_score"]["param_value"]
                else:
                    bitscore = ''

                nucl_seq = entries[entry]["dna_sequence"]["sequence"]
                seqlen = len(nucl_seq)
                organism = entries[entry]["NCBI_taxonomy"]["NCBI_taxonomy_name"]
                product = json_db[model]["ARO_description"]
                ref_db = "CARD-{}".format(db_version)

                append_json(meta_data, gene_id, db=ref_db, bitscore=bitscore, \
                            genefam=acc, product=product, organism=organism, \
                            seqlen=seqlen, gene=aro_name)

            # Generate CARD fasta files by model type, if requested
            if args.card_fasta:
                model_type = json_db[model]["model_type"]

                try:
                    out_fa = out_map[model_type]["outfile"]
                except KeyError:
                    continue
    
                seq_type = "protein_sequence" if "protein" in model_type else \
                           "dna_sequence"

                sequence = entries[entry][seq_type]["sequence"]
                out_fa.write(">{}\n{}\n".format(gene_id, sequence))

    # Output internal ontology
    for ontology in ontologies:
        out_h.write("{}\t{}\n".format(ontology, reaction_alts(ontologies[ontology])))

    # Output internal gene metadata
    if args.card_genes:
        json.dump(meta_data, args.card_genes, sort_keys=True, indent=4, separators=(',', ': '))

    # Database statistics:
    print("Total database size: {}".format(db_totals), file=sys.stderr)
    for model_type in out_map:
        print("    models in '{}':\t{!s}".format(model_type, out_map[model_type]["counts"]), file=sys.stderr)


def sub_kegg(args):
    """
    Subcommand for generating internal reference files for the KEGG database
    """
    pass


def sub_uniprot(args):
    """
    Subcommand for generating internal reference files for the uniProt database
    """
    pass


_hierarchy = {
    "antibiotic resistant gene variant or mutant": 1, 
    "gene conferring resistance via absence": 1, 
    "gene involved in antibiotic sequestration": 1, 
    "antibiotic target protection protein": 1, 
    "antibiotic target modifying enzyme": 1.5, 
    "antibiotic resistance gene cluster, cassette, or operon": 1, 
    "antibiotic inactivation enzyme": 1.5, 
    "protein modulating permeability to antibiotic": 2, 
    "protein(s) conferring antibiotic resistance via molecular bypass": 2, 
    "antibiotic target replacement protein": 2, 
    "gene involved in self-resistance to antibiotic": 2, 
    "gene altering cell wall charge": 2, 
    "protein(s) and two-component regulatory system modulating antibiotic efflux": 3, 
    "gene modulating beta-lactam resistance": 3, 
    "efflux pump complex or subunit conferring antibiotic resistance": 4
             }


def sort_by_level(ontology, hierarchy=_hierarchy):
    """
    Sort CARD categories in ascending order using a custom hierarchy
    """
    res_re = re.compile("(?<=determinant of ).+(?= resistance)")
    anti_re = re.compile("(?<=determinant of resistance to ).+(?= antibiotics)")

    levels = []
    determs = []
    others = []
    # Obtain list of category levels
    for component in ontology:
        try:
            levels.append(hierarchy[component])
            others.append(component)
        except KeyError:
            if "determinant of" in component:
                determs.append(component)
            else:
                print("Category '{}' does not yet exist in custom hierarchy. Please "
                      "rectify by adding to the hierarchy before attempting to "
                      "generate the internal ontology".format(';'.join(ontology)), file=sys.stderr)
                sys.exit(1)

    # Sort categories in ascending order by level
    order = []
    for level in sorted(levels):
        order.append(levels.index(level))

    ontology = [others[i] for i in order]

    # Check if "determinant of ..." is one of the categories
    if determs:
        determs = sorted(determs)
        # Does determinant confer multi-drug resistance?
        if len(determs) > 1:
            try:
                cat_start = [anti_re.search(i).group() for i in determs[:-1]]
                cat_end = anti_re.search(determs[-1]).group()
            except AttributeError:
                cat_start = [res_re.search(i).group() for i in determs[:-1]]
                cat_end = res_re.search(determs[-1]).group()

            determ = "determinant of multi-drug resistance to {} and {} "\
                       "antibiotics".format(', '.join(cat_start), cat_end)

        else:
            determ = determs[0]

        ontology.append(determ)

    return ontology


def append_json(json_data, acc, db='', bitscore='', gene='', genefam='', organism='', product='', seqlen=''):
    """
    Append an entry to a JSON database of gene metadata
    """
    json_data[acc] = {
                      'organsim': organism,
                      'product': product,
                      'gene_length': seqlen,
                      'gene': gene,
                      'gene_family': genefam,
                      'database': db,
                      'bitscore_threshold': bitscore
                     }

    return(json_data)


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


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parent_parser = argparse.ArgumentParser(add_help=False)
    subparsers = parser.add_subparsers(title="subcommands",
        help="exactly one of these commands is "
             "required")
    # FOAM-specific arguments
    foam_parser = subparsers.add_parser('foam',
        parents=[parent_parser],
        help="generate internal files for the FOAM ontology")
    foam_parser.add_argument('-f', '--foam-in',
        metavar='INFILE',
        dest='foam_in',
        action=Open,
        mode='rb',
        help="input tabular FOAM ontology file relating Kegg Orthologs (KOs) "
             "to a hierarchical representation of biochemical functions and "
             "pathways")
    foam_parser.add_argument('-fo', '--foam-out',
        metavar='OUTFILE',
        dest='foam_out',
        action=Open,
        mode='wt',
        default=sys.stdout,
        help="use the FOAM ontology to output a tabular KO to pathway mapping "
             "file, formatted for internal use [default: output to stdout]")
    foam_parser.set_defaults(func=sub_foam)
    # CARD-specific arguments
    card_parser = subparsers.add_parser('card',
        parents=[parent_parser],
        help="generate internal files for the CARD reference database. CARD is not as established a database as KEGG or UniProt. Future updates may break compatibility, requiring modifications to the program")
    card_parser.add_argument('-c', '--card-json',
        metavar='INFILE',
        dest='card_in',
        action=Open,
        mode='rb',
        help="input complete CARD database in JSON format")
    card_parser.add_argument('-co', '--card-out',
        metavar='OUTFILE',
        dest='card_out',
        action=Open,
        mode='wt',
        default=sys.stdout,
        help="use the Antibiotic Resistance Ontology (ARO) to output a tabular "
             "AR category mapping file, formatted for internal use [default: "
             "output to stdout]")
    card_parser.add_argument('-cg', '--card-genes',
        metavar='OUTFILE',
        dest='card_genes',
        action=Open,
        mode='wt',
        help="output gene metadata file in JSON format. Metadata file "
             "contains accessions along with their corresponding gene/gene "
             "family ids, gene lengths, product and organismal information, "
             "etc.")
    card_parser.add_argument('-cs', '--card-snp',
        metavar='OUTFILE',
        dest='card_snp',
        action=Open,
        mode='wt',
        help="output tabular file ")
    card_parser.add_argument('-cm', '--card-models',
        metavar='MODEL [,MODEL]',
        dest='card_models',
        default=["protein homolog model", "protein variant model", "rRNA gene variant model", "protein overexpression model"],
        help="Output information on the following CARD model types [default: "
              "'protein homolog model', 'protein variant model', 'rRNA gene "
              "variant model', 'protein overexpression model']")
    card_parser.add_argument('--card-fastas',
        dest='card_fasta',
        action='store_true',
        help="")
    card_parser.set_defaults(func=sub_card)
    # KEGG-specific arguments
    kegg_parser = subparsers.add_parser('kegg',
        parents=[parent_parser],
        help="generate internal files for KEGG. Currently not-supported")
    kegg_parser.add_argument('-k', '--kegg-in',
        metavar='INFILE',
        dest='kegg_in',
        action=Open,
        mode='rb',
        help="input DBGET flat file containing KO entries")
    kegg_parser.add_argument('-km', '--kegg-modules',
        metavar='INFILE',
        dest='kegg_modules',
        action=Open,
        mode='rb',
        help="input DBGET flat file containing MODULE entries")
    kegg_parser.add_argument('-ko', '--kegg-out',
        metavar='OUTFILE',
        dest='kegg_out',
        action=Open,
        mode='wt',
        default=sys.stdout,
        help="use the KEGG ontology (ko) to output a tabular KO to pathway "
             "mapping file, formatted for internal use [default: output to "
             "stdout]")
    kegg_parser.add_argument('-kg', '--kegg-genes',
        metavar='OUTFILE',
        dest='kegg_genes',
        action=Open,
        mode='wt',
        help="output gene metadata file in JSON format. Metadata file "
             "contains accessions along with their corresponding gene/gene "
             "family ids, gene lengths, product and organismal information, "
             "etc.")
    kegg_parser.set_defaults(func=sub_kegg)
    # UniProt-specific arguments
    uniprot_parser = subparsers.add_parser('uniprot',
        parents=[parent_parser],
        help="generate internal files for the UniProt protein databases")
    uniprot_parser.add_argument('-u', '--uniprot-in',
        metavar='INFILE',
        dest='uniprot_in',
        action=Open,
        mode='rb',
        help="input DBGET flat file containing UniProt entries")
    uniprot_parser.set_defaults(func=sub_uniprot)
    args = parser.parse_args()

    args.func(args)


if __name__ == "__main__":
    main()
    sys.exit(0)
