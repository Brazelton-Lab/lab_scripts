#! /usr/bin/env python
# converts OTU file from QIIME into .list file for mothur
# incorporates sequence names from provided .names file

# usage: python convert_qiime_otus_to_mothur_list.py otus.txt VAMPS_JGIbtrim30.unique.good.filter.names 31776
# user must provide the following input manually:

dist = '0.03'

import sys
infilename = sys.argv[1]
namesfilename = sys.argv[2]
otus = sys.argv[3]
outfilename = infilename + '.list'
outfile = open(outfilename, 'a')
outfile.write(dist+'\t'+otus+'\t')

namesfile = open(namesfilename)
N=[]
for nline in namesfile:
        nline = nline.strip('\n')
        ntabs = nline.split('\t')
        names = ntabs[1]
        #names = names.split(',')
        N.append(names) 

infile = open(infilename)
for line in infile:
        tabs = line.split('\t')
        seq = tabs[1].strip('\n')
        #print seq
        found = 'no'
        for names in N:
                #print names
                if seq in names: 
                        outfile.write(names)
                        found = 'yes'
        if found == 'no': 
                print seq,
                print 'not found in names file'         
        if len(tabs) >2:
                #print tabs
                found = 'no'
                for tab in tabs[2:]:
                        outfile.write(',')
                        seq=tab.strip('\n')
                        for names in N:
                                if seq in names: 
                                        outfile.write(names)
                                        found = 'yes'
                                        #print names
                        if found == 'no': 
                                print seq,
                                print 'not found in names file'                 
        outfile.write('\t')
infile.close()
outfile.close()
namesfile.close()