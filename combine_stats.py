#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Copyright:

    combine_stats  Combines all the statistical data from multiple sources for each assembly into one file.

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
__version__ = '1.1'


import pandas as pd
import re
import argparse


def main():
    ap = argparse.ArgumentParser(description="Combines all the statistical data from multiple sources for each assembly into one file.")
    ap.add_argument("-a", dest="ale", nargs="+", required=True, help="One or more input .ale files with or without the header info.")
    ap.add_argument("-b", dest="bowtie2", nargs="+", required=True, help="One or more input bowtie2 files.")
    ap.add_argument("-m", dest="metaquast", nargs="+", required=False, help="One or more input metaquast files.")
    ap.add_argument("--other", nargs="+", required=False, help="One or more input csv files with the same name index.")
    ap.add_argument("-o", dest="output", default="combined_stats.csv", help="The file output name.")
    args = ap.parse_args()

    combine_stats(args.output, args.ale, args.bowtie2, args.metaquast, other=args.other)

    
def combine_stats(output, ale_file, bt2_file, metaquast_file=None, other=None,
                  a_id="file_name", b_id="file_name", m_id="Assembly"):
    ale = merge_files(ale_file, a_id)
    bt2 = merge_files(bt2_file, b_id)
    
    ale_re = re.compile("([\.\w]+\/)*(\w+)((_ale_data)|(\.ale))")
    bt2_re = re.compile("([\.\w]+\/)*(\w+)[_\/]?mapping\.log")
    create_id(ale, ale_re, a_id)
    create_id(bt2, bt2_re, b_id)
    
    mqt = None
    if metaquast_file:
        data = pd.read_csv(metaquast_file[0], sep="\t")  # the first file
        for file in metaquast_file[1:]: # the rest of the files
            temp = pd.read_csv(file, sep="\t")
            data = data.append(temp)
        mqt = pd.DataFrame(data)
        mqt_re = re.compile("(.*)")  # everything
        create_id(mqt, mqt_re, m_id, reg_group=1)
    
    results = ale.merge(bt2, on="id")
    if metaquast_file:
        results = results.merge(mqt, on="id")
    if other:
        for o in other:
            extra = pd.read_csv(o)
            results = results.merge(extra, on="id")
    results.to_csv(output, index=False)
    
    
def create_id(data, reg, previous_index, reg_group=2):
    """Extracts the id from the previous_index column using a regex then sets the new id as the index."""
    names = list()
    for a in data[previous_index]:
        am = reg.match(a)
        if not am:
            print("Warning: skipping invalid file name: " + str(a))
            continue
        names.append(am.group(reg_group))
    data["id"] = names
    data.set_index("id")
    data.drop(previous_index, axis=1, inplace=True)
    
    
def merge_files(files, index_col=None, sep=","):
    """Combines multiple csv's with similar data into one organized by the index_col."""
    data = pd.read_csv(files[0], sep=sep)  # the first file
    if index_col:
        for file in files[1:]: # the rest of the files
            temp = pd.read_csv(file, sep=sep)
            data = data.merge(temp, on=index_col, how="outer")
    data = pd.DataFrame(data)
    return data
    
    
if __name__ == "__main__":
    main()
    
