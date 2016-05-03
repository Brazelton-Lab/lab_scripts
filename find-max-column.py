#! /usr/bin/env python

# report column with the highest number
# assumes comma-separated-values and unix linebreaks
# usage:
# python find-max-column.py final.an.0.03.subsample.shared.nosingles.transposed.csv

import sys
filename = sys.argv[1]
file = open(filename)
outfilename = filename.replace('.csv','.maxcolumn.csv')
outfile = open(outfilename, 'a')

for line in file:
        line = line.strip('\n')
        headers = line.split(',')
        headers = headers[1:]
        break
        
for line in file:
        line = line.strip('\n')
        columns = line.split(',')
        OTU = columns[0]
        numbers = []
        for i in columns[1:]: numbers.append(float(i))
        maxcolumn = numbers.index(max(numbers))
        outfile.write(OTU)      
        outfile.write(',')
        outfile.write(headers[maxcolumn])
        outfile.write(',')
        outfile.write(str(max(numbers)))
        outfile.write('\n')
        
file.close()
outfile.close() 
