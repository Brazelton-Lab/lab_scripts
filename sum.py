#! /usr/bin/env python3

# simple python script to sum values in each column because I was scared of awk not giving me the values in the proper order
# assumes comma separated values


import sys
import csv
import pandas as pd

infile = sys.argv[1]
outfile = infile +'.sums.csv'

l = []
with open(infile) as f:	
	# read data as a dataframe with pandas	
	df = pd.read_csv(f, index_col=[0])
	total = df.sum()

total.to_csv(outfile)