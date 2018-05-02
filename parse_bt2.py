#!/usr/bin/env python3.6
# -*- coding: utf-8 -*-
"""

Copyright:

    parse_bt2  Extracts the statistics found in bowtie2 output files to a csv.

    Copyright (C) 2016  William Brazelton, Nickolas Lee

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


__author__ = 'Nickolas Lee'
__license__ = 'GPLv3'
__status__ = 'Production'
__version__ = '1.0'


import re
import pandas as pd
import argparse


def main():
    ap = argparse.ArgumentParser(description="Extracts the statistics found in bowtie2 output files to a csv.")
    ap.add_argument("-o", "--output_file", default="bt2.csv", help="The name of the csv file to be written.")
    ap.add_argument("files", nargs="+", help="A list of bowtie2 output files.")
    args = ap.parse_args()
    
    data = pd.DataFrame()
    for file in args.files:
        bt = parse_bt(file)
        for k in bt.keys():
            data.loc[file,k] = bt[k]
    print("Saving to: "+ args.output_file)
    data.to_csv(args.output_file, index=False)
    
    
def parse_bt(file):
    """Extracts the statistics found in a bowtie2 output file."""
    nums = list()
    start_re = re.compile("\d+ reads; of these:")
    amount_re = re.compile("^\s*(\d+)\s")
    perc_re = re.compile("(\d+\.\d+)\%")
    with open(file) as bt:
        start = False
        al_rate = 0
        for line in bt.readlines():
            if start_re.match(line):
                start = True
            if start:
                amount = amount_re.match(line)
                percent = perc_re.match(line)
                if percent:
                    al_rate = percent.group(1)
                if amount:
                    nums.append(amount.group(1))
    col_names = ["Total reads", "Total paired reads", "Paired reads with no concordant alignment",
                 "Uniquely aligned paired reads", "Paired reads with multiple alignments",
                 "", "Uniquely discordantly aligned paired reads", "Paired reads with no alignment", 
                 "","","","", "Total single reads", "Unaligned single reads", 
                 "Uniquely aligned single reads", "Single reads with multiple alignments"]
    data = dict()
    if nums:
        for n in range(len(nums)):
            data[col_names[n]] = nums[n]
        data.pop("")  # removes unneeded key
        data["Total mapped fragments"] = int(data["Total paired reads"]) - int(data["Paired reads with no alignment"])
        if float(al_rate) > 0:
            data["Overall alignment rate"] = al_rate
    data["file_name"] = file
    return data

    
if __name__ == "__main__":
    main()
