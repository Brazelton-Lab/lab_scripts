#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Copyright:

    parse_ale  Parses .ale files and produces mapped contig coverage statistics for each assembly. 

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
__maintainer__ = __author__
__status__ = 'Production'
__version__ = '2.0'


import pandas as pd
import re
import matplotlib.pyplot as pp
from statistics import mean, stdev, mode, median, StatisticsError
from scipy.stats import skew
import argparse
import pickle


def main():
    ap = argparse.ArgumentParser(description="Parses .ale files and produces mapped contig coverage statistics for each assembly. To save compute time on repeated analyses, it can also produce and consume intermediate .pkl files made from the assembly data.")
    ap.add_argument("-o", "--output_stats", default="ale_statistics.csv", help="Name of the output csv.")
    ap.add_argument("-p", "--plot", action="store_true", default=False, help="If specified then save histogram plots of data for each input file to the working directory.")
    ap.add_argument("--header_only", action="store_false", default=True, help="If specified then do not compute coverage statistics, plots, or intermediate files.")
    ap.add_argument("-s", "--save_data", action="store_true", default=False, help="Whether to save intermediate data.")
    ap.add_argument("-i", "--intermediates", action="store_true", default=False, help="Specifies that the input files are intermediate .pkl data files not .ale files.")
    ap.add_argument("-v", "--verbose", action="store_true", default=False, help="Print out which step is being worked on. Errors are always printed.")
    ap.add_argument("files", nargs="+", help="A list of ale files.")
    args = ap.parse_args()
    analyse_ale(args.files, plot=args.plot, coverage=args.header_only, out_file_name=args.output_stats, save=args.save_data, is_pickle=args.intermediates, verbose=args.verbose)
    
    
def analyse_ale(files, plot = False, coverage = True, out_file_name = "ale_statistics.csv", save=False, is_pickle=False, verbose=True):
    """Parses .ale files and produces mapped contig coverage statistics for each assembly. To save compute time on repeated analyses, it can also produce and consume intermediate .pkl files made from the assembly data."""
    data = pd.DataFrame()
    c = -1
    for file in files:
        c += 1
        if not is_pickle:
            file_re = re.compile("(\.\.\/)*(data\/)?(singles\/)?(\w+\/)*(\w+).ale")
            m = file_re.match(file)
            file_id = m.group(5)
            if verbose:
                print(file_id)
            if file_id == "":
                print("Warning! Skipping invalid file name: "+file)
                continue
            data.loc[c,"file_name"] = file_id
            gath,header_end = gather_stats(file)
            for g in gath.items():
                data.loc[c,g[0]] = g[1]
        else:
            file_re = re.compile("(\w+\/)*(\w+)(_ale_data)?\.pkl")
            m = file_re.match(file)
            file_id = m.group(2)
            if file_id == "":
                print("Warning! Skipping invalid file name: "+file)
                continue
            data.loc[c,"file_name"] = file_id
            
        if coverage:
            if verbose:
                print("Reading coverage data")
                
            if not is_pickle:
                contig_coverages = read_coverage(file, header_end)
                
                if save:
                    j_name = file_id+"_ale_data.pkl"
                    if verbose:
                        print("Saving to "+j_name)
                    with open(j_name, "wb") as output:
                        pickle.dump(contig_coverages, output)
            else:
                with open(file, "rb") as put:
                    contig_coverages = pickle.load(put)
            
            if verbose:
                print("Computing statistics")
            
            avgs = dict()
            maxs = dict()
            for a in contig_coverages.items():
                avgs[a[0]] = dict_avg(a[1])
                maxs[a[0]] = dict_max(a[1])
            data.loc[c,"Avg contig avg bp coverage"] = mean(avgs.values())
            data.loc[c,"Avg contig max bp coverage"] = mean(maxs.values())
            
            total_coverage, lengths = compute_total_coverage(contig_coverages)

            data.loc[c,"Min contig total coverage"] = min(total_coverage)
            data.loc[c,"Avg contig total coverage"] = mean(total_coverage)
            data.loc[c,"Stdev of contig total coverage"] = stdev(total_coverage)
            data.loc[c,"Max contig total coverage"] = max(total_coverage)
            try:
                data.loc[c,"Mode of contig total coverage"] = mode(total_coverage)
            except StatisticsError as se:
                if verbose:
                    print("Error computing mode: "+str(se))
                pass
            total_coverage.sort()
            data.loc[c,"Median of contig total coverage"] = median(total_coverage)
            data.loc[c,"Skew of contig total coverage"] = skew(total_coverage, nan_policy="omit")
            
            data.loc[c,"Min contig length"] = min(lengths)
            data.loc[c,"Avg contig length"] = mean(lengths)
            data.loc[c,"Stdev of contig lengths"] = stdev(lengths)
            data.loc[c,"Max contig length"] = max(lengths)
            try:
                data.loc[c,"Mode of contig length"] = mode(lengths)
            except StatisticsError as se:
                if verbose:
                    print("Error computing mode: "+str(se))
                pass
            lengths.sort()
            data.loc[c,"Median of contig length"] = median(lengths)
            data.loc[c,"Skew of contig length"] = skew(lengths, nan_policy="omit")
            
            assembly = summarize_assembly(contig_coverages)

            data.loc[c,"Assembly min bp coverage"] = dict_min(assembly)
            data.loc[c,"Assembly avg bp coverage"] = dict_avg(assembly)
            data.loc[c,"Assembly max bp coverage"] = dict_max(assembly)
            
            if plot:
                pp.title("Histogram of mapping coverage by amount of basepairs for "+file_id)
                pp.xlabel("Coverage")
                pp.ylabel("Basepair count")
                pp.scatter([cov for cov in assembly.keys()], [bp for bp in assembly.values()])
                histo = file_id+"_plot.png"
                if verbose:
                    print("Saving histogram to " + histo)
                pp.savefig(histo)
                pp.close()
    
    if data.shape[0] != 0:  # if not empty
        try:
            data.drop("Reference", axis=1, inplace=True)  # I don't know what this column is for
        except ValueError:
            pass
        if verbose:
            print("Saving data to " + out_file_name)
        data.to_csv(out_file_name, index=False)
    else:
        print("Error: Nothing to report")

        
def gather_stats(file):
    """Parses the header of an .ale file."""
    stats = dict()
    reg = re.compile("\# (.+): (.+)\n?")
    with open(file) as fi:
        count = 0
        for line in fi:
            parts = reg.match(line)
            if parts:
                    stats[parts[1]] = parts[2]
                    count += 1
            else:
                break
    return stats, count


def read_coverage(file, start_line):
    """Produces a dictionary containing mapped coverage per nucleotide as keys with a count 
    of their presence as values (a histogram dictionary)."""
    with open(file) as fi:
        run = dict()
        for st in range(start_line + 1):
            fi.readline()
        for line in fi:
            line = line.replace("\n","")
            if len(line) == 0:
                continue
            nums = line.split(" ")
            if len(nums) != 7:  # invalid input
                continue
            val = int(nums[2])
            try:
                run[nums[0]]
            except KeyError:
                run[nums[0]] = dict()
            try:
                run[nums[0]][val] += 1  # count
            except KeyError:
                run[nums[0]][val] = 1  # initialize value
        return run
            
    
def dict_avg(di):
    """Finds the unique item average of a histogram dictionary"""
    if len(di) > 0:  # avoids divide by zero
        result = 0
        denom = 0
        for d in di.items():
            result += int(d[0]) * int(d[1])
            denom += int(d[1])  # the total number of bp
        return float(result) / denom


def dict_min(di):
    """Finds the smallest key out of a histogram dictionary"""
    begining = True
    result = None
    for v in di.items():
        if begining:
            result = v[0]
            begining = False
        if v[0] < result:
            result = v[0]
    return result
        

def dict_max(di):
    """Finds the largest key out of a histogram dictionary"""
    begining = True
    result = None
    for v in di.items():
        if begining:
            result = v[0]
            begining = False
        if v[0] > result:
            result = v[0]
    return result
    

def compute_total_coverage(di):
    """Compute the sum of the coverage for every bp and the length for every contig."""
    total_coverage = list()
    lengths = list()
    for contig in di.values():
        length = 0
        area = 0
        for cc in contig.items():
            area += cc[0] * cc[1]   # sum(coverage * count) = total coverage area 
            length += cc[1]         # count of nucleotides
        total_coverage.append(area)
        lengths.append(length)
    return total_coverage, lengths


def summarize_assembly(di):
    """Add all the histogram dictionaries into one"""
    assembly = dict()
    for cont in di.values():
        for a in cont.items():
            if len(a) > 0:
                b = int(a[0])
                try:
                    assembly[b]
                except KeyError:
                    assembly[b] = 0
                assembly[b] += a[1]
    return assembly
                
    
if __name__ == "__main__":
    main()
