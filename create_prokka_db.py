#! /usr/bin/env python
"""
create a database fasta for use with prokka using m5nr fasta and translation 
files
"""

import sys
import os
import argparse
import re

def io_check(infile, mode='rU'):
    try:
        fh = open(infile)
    except IOError as e:
        try:
            fh = open(infile, mode)
        except IOError as e:
            print(e)
            sys.exit(1)
    else:
        if mode == 'w':
            print("{} already exists in cwd".format(os.path.basename(infile)))
            sys.exit(1)
        fh.close()
    return infile

def parse_fasta(fasta):
    sequences = {}
    with open(fasta) as in_h:
        line = in_h.readline()
        while line:
            line = line.strip()
            if line.startswith('>'):
                m5nr = line.strip('>\n')
                m5nr = m5nr.split()[0]
            else:
                try:
                    sequences[m5nr]['sequence'] += line
                except KeyError:
                    sequences[m5nr] = {'sequence': line}
            line = in_h.readline()
    return sequences

def main():
    out_db = io_check(args.out, 'w')
    out_mapper = io_check("function_idmapper.csv", 'w')
    seqs = parse_fasta(args.fasta)
    genes = {}

    func_map = args.func
    with open(func_map) as in_h:
        r = re.compile("((?<=\()EC.+(?=\)))")
        for line in in_h:
            line = line.strip().split('\t')
            m5nr = line[0]
            gene = line[1]
            taxon = line[3]
            matched = r.search(line[2])
            if not matched:
                product = line[2]
                ec = ''
            else:
                ec = matched.group()
                product = line[2].replace('({})'.format(ec), '').rstrip()
            seqs[m5nr]['gene'] = gene
            seqs[m5nr]['ec'] = ec
            seqs[m5nr]['product'] = product
            seqs[m5nr]['taxon'] = taxon
            genes[gene] = m5nr
    
    ortholog_map = args.go
    with open(ortholog_map) as in_h:
        for line in in_h:
            line = line.strip().split('\t')
            m5nr = line[0]
            orth = line[1]
            try:
                seqs[m5nr]['gene ortholog'] = orth
            except KeyError:
                print("no match for ortholog: {}".format(orth))
                sys.exit(1)

    with open(out_db, 'w') as db_h:
        with open(out_mapper, 'w') as mapper_h:
            for gene in sorted(genes):
                seq_id = genes[gene]
                ec = seqs[seq_id]['ec']
                product = seqs[seq_id]['product']
                try:
                    go = seqs[seq_id]['gene ortholog']
                except KeyError:
                    continue
                seq = seqs[seq_id]['sequence']
                output = ">{} {}~~~{}~~~{}\n{}\n".format(go, ec, gene, product, seq)
                db_h.write(output)

                seq_len = 3 * len(seq)
                taxon = seqs[seq_id]['taxon']
                output = "{}\t{}\t{}\t{}\n".format(gene, go, str(seq_len), taxon)
                mapper_h.write(output)

if __name__ == "__main__":
   parser = argparse.ArgumentParser(description="create a prokka database")
   parser.add_argument('fasta', metavar='<fasta db>',
                       type=io_check,
                       help="reference database in fasta format")
   parser.add_argument('go', metavar='<gene families>',
                       type=io_check,
                       help="m5nr id to gene orthology mapper file")
   parser. add_argument('func', metavar='<genes>',
                        type=io_check,
                        help="m5nr id to gene mapper file")
   parser.add_argument('-o', '--out', metavar='FILE',
                       default='db',
                       help="output file name")
   args = parser.parse_args()
   main()
   sys.exit(0)
