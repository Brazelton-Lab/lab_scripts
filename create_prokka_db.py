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
                m5nr = line.strip('>')
                m5nr = m5nr.split()[0]
            else:
                seq_len = len(line)
                try:
                    sequences[m5nr]['sequence'] = sequences[m5nr].get('sequence', '') + line
                    sequences[m5nr]['length'] = sequences[m5nr].get('length', 0) + seq_len
                except KeyError:
                    sequences[m5nr] = {'sequence': line, 'length': seq_len}
            line = in_h.readline()
    return sequences

def main():
    out_db = io_check(args.out)
    out_mapper = io_check("function_idmapper.csv")
    seqs = parse_fasta(args.fasta)

    ortholog_map = args.go
    with open(ortholog_map) as in_h:
        for line in in_h:
            line = line.strip().split('\t')
            m5nr = line[0]
            orth_ident = line[1]
            try:
                seqs[m5nr]['gene ortholog'] = orth_ident
            except KeyError:
                continue

    func_map = args.func
    with open(func_map) as in_h:
        r = re.compile("((?<=\()EC.+(?=\)))")
        for line in in_h:
            line = line.strip().split('\t')
            m5nr = line[0]
            gene_ident = line[1]
            taxon = line[3]
            matched = r.search(line[2])
            if not matched:
                product = line[2]
                ec_num = ''
            else:
                ec_num = matched.group()
                product = line[2].replace('({})'.format(ec_num), '').rstrip()
            try:
                seqs[m5nr]['gene'] = gene_ident
                seqs[m5nr]['ec'] = ec_num
                seqs[m5nr]['product'] = product
                seqs[m5nr]['taxon'] = taxon
            except KeyError:
                print("mismatch: {}".format(ident))

    with open(out_db, 'w') as db_h:
        with open(out_mapper, 'w') as mapper_h:
            for seq_id in seqs:
                try:
                    gene = seqs[seq_id]['gene']
                    ec = seqs[seq_id]['ec']
                    product = seqs[seq_id]['product']
                except KeyError:
                    continue
                seq = seqs[seq_id]['sequence']
                seq_len = seqs[seq_id]['length']
                output = ">{} {}~~~{}~~~{}\n{}\n".format(gene, ec, seq_id, product, seq)
                db_h.write(output)
                try:
                    go = seqs[seq_id]['gene ortholog']
                    taxon = seqs[seq_id]['taxon']
                except KeyError:
                    continue
                seq_len = seqs[seq_id]['length']
                output = "{}\t{}\t{}\t{}\n".format(gene, go, seq_len, taxon)
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
