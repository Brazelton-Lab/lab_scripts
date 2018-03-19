#!/usr/bin python3
# -*- coding: utf-8 -*-
"""
Created on 2018-03-15

@author: Nickolas Lee

Copyright:

    compute_usage  Downloads the Slurm history database and produces a csv report and/or plots showing program useage statistics. 
    Copyright (C) 2016  Nickolas Lee, William Brazelton

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
__version__ = '1.1.0'


import pandas as pd
import matplotlib.pyplot as pp
import re
import scipy.stats
import argparse


def main():
    parser = argparse.ArgumentParser(description="Reads "+
            "the delimited output of slurm and produces a "+
            "csv report and/or plots of past program useage statistics. "+
            "First obtain the slurm CSV by running this command on the cluster: "+
            "sacct -ap -o ALL --noconvert -S \"2015-01-01\" > slurm_output.csv")
    parser.add_argument("-v", "--version", action='version', version=__version__)
    parser.add_argument("-c", "--clean", action='store_true', default="false",
                        help="Cleans the data by removing all rows containing NAN values.")
    parser.add_argument("-r", "--make_report", action='store_true', default="false",
                        help="Writes statistics to a CSV file.")
    parser.add_argument("-p", "--make_plot", action='store_true', default="false",
                        help="Creates plots of the data for more refined analysis.")
    parser.add_argument("-o", "--output_folder", action='store', default="",
                        help="Specify the output folder. Must end with a slash character.")
    parser.add_argument("-q", "--query", action='store', default="",
                        help="Specify the name of the program to look for. The default is all programs.")
    parser.add_argument("-d", "--data", action='store', default="slurm_output.csv",
                        help="Give the location of the slurm CSV output.")
    args = parser.parse_args()
    
    compute_usage(program_name = args.query, folder = args.folder, 
                   data_file = args.data,
                   clean=args.clean, make_report=args.make_report, 
                   make_plot=args.make_plot)
    
    
def compute_usage(program_name = "", folder = "", data_file = "",
                   clean = False, make_report = False, make_plot = False):
    """Reads in the delimited output of slurm and produces a csv report and/or 
    plots of useage statistics. One program can be specified using program_name. 
    But, the default uses all programs. Optionally the data can be cleaned 
    before use by removing all rows with NAN values."""
    
    col = ["JobName", "MaxRSS", "AllocCPUS", "Elapsed", "State", "MaxDiskRead",
           "MaxDiskWrite"]
    
    file = pd.read_csv(data_file, usecols=col, delimiter="|")
    file = pd.DataFrame(file)
    if clean:
        file.dropna(inplace=True)  # hopefully removes most of the incomplete runs
        clean_name = "cleaned_"
    else:
        file.fillna(0, inplace=True)
        clean_name = ""
    
    # clean up the formatting
    file["MaxRSS"] = file["MaxRSS"].apply(SI_Unit_convert)
    file["MaxDiskRead"] = file["MaxDiskRead"].apply(SI_Unit_convert)
    file["MaxDiskWrite"] = file["MaxDiskWrite"].apply(SI_Unit_convert)
    file["Elapsed"] = file["Elapsed"].apply(fix_time_format)
    file["Elapsed"] = pd.to_timedelta(file["Elapsed"]) 
    
    # run query
    if program_name:
        file = file[file["JobName"] == program_name]
    file = file[file["State"] == "COMPLETED"]
    if file.shape[0] == 0: 
        print("Error: No data fit these conditions")
        return

    if program_name:
        program_name = program_name+"_"
        
    # produce CSV report
    if make_report:
        jobs = file.groupby(by="JobName")
        des = jobs.describe()
        des.to_csv(folder+program_name+clean_name+"completed_job_statistics.csv")

    # product plots
    if make_plot:
        # remove columns that do not contain numbers
        ana = [c for c in col if c != "JobName"]
        ana = [c for c in ana if c != "State"]
        
        xaxis = ["MaxDiskRead", "AllocCPUS"]
        for xl in xaxis:
            for by in ana:
                if xl != by:
                    pp.scatter(file[xl], pd.to_numeric(file[by]))
                    pp.xlabel(xl)
                    pp.ylabel(by)
                    # calculate the best fit line of the points for extrapolation
                    try:
                        slope, intercept, r_value, pvalue, sterr = scipy.stats.linregress(file[xl], 
                                                          pd.to_numeric(file[by]))
                        equ = "y=("+str(format(slope, ".3"))+")x+"+str(int(intercept))
                        pp.plot(file[xl], slope*file[xl] + intercept, "r")
                        pp.legend([equ, "R="+str(format(r_value, ".3"))])
                    except (ValueError):  # this is for nan values
                        pass
                    pp.savefig(folder+program_name+clean_name+"completed_"+xl
                                       +"_by_"+by+".png")
                    pp.close()  # prevents plots from forming on top of other plots


def SI_Unit_convert(data):
    """Converts strings with SI unit notation into numbers."""
    data = str(data)
    if data and data != "nan":
        num = 0
        try:
            num = float(data[:-1])
            if data[-1] == "K":
                num = num*1000
            if data[-1] == "M":
                num = num*1000000
            if data[-1] == "G":
                num = num*1000000000
        except (ValueError):
            pass
        return num
    
    
def fix_time_format(data):
    """Converts days in the Slurm output into hours for pandas."""
    day_re = re.compile("(\d+)\-(\d\d)(\:\d\d\:\d\d)")
    day = day_re.match(data)
    if day:
        hours = int(day.group(1)) * 24 + int(day.group(2))
        if hours >= 0:  # for the case of integer overflow and negative deltas
            return str(hours) + day.group(3)
    return fix_PM(data)


def fix_PM(data):
    """Adds 12 hours to a time ending in PM."""
    pm_re = re.compile("(\d\d)(\:\d\d\:\d\d) PM")
    pm = pm_re.match(data)
    if pm:
        hours = int(pm.group(1)) + 12
        return str(hours) + pm.group(2)
    return data


if __name__ == "__main__":
    main()
