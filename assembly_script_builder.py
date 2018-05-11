#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Copyright:

    assembly_script_builder  Creates shell scripts that perform assembly of the samples and computes statistics for compareing assemblies.
    
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
__version__ = '1.3'


import argparse
import re


def main():
    ap = argparse.ArgumentParser(description='Creates shell scripts that perform assembly of the samples and computes statistics for compareing assemblies. By default it automatically groups files by filename before the first period. Assumes the input files have already been preprocessed and labled with: sample name, period, type of read file, period, and then other file endings.')
    ap.add_argument("-s", '--samples', nargs="+", required=True, help='list of preprocessed read input files')
    ap.add_argument("-g", "--group_ids", nargs="+", help="a list of ordered group ids corresponding to how the input files should be grouped together for assembly")
    ap.add_argument("-d", "--data_folder", default=".", help="the output folder for data")
    ap.add_argument("-r", "--results_folder", default=".", help="the output folder for results")
    ap.add_argument("-o", "--script_folder", default=".", help="the output folder for scripts")
    ap.add_argument("-l", "--log_folder", default=".", help="the output folder for logs during runtime")
    ap.add_argument("-c", "--cpus", type=int, default=32, help="the number of cpus to use")
#    ap.add_argument("-a", "--assembler", default="ray", help="the assembler to use for assembly. The default is Ray meta designated 'ray'. This is currently the only option.")
    ap.add_argument("-p", "--plot", action="store_true", default=False, help="whether to save plots from the analysis")
    ap.add_argument("--pool", action="store_true", default=False, help="assemble all files in one assembly together")
    ap.add_argument("--separate", action="store_true", default=False, help="assemble all files in different assemblies. Make sure each file name has a unique id before the first period otherwise the script files generated will be overwritten.")
    args = ap.parse_args()
    
    # prepare variables    
    if args.script_folder and args.script_folder[-1] != "/":
        args.script_folder = args.script_folder + "/"
    
    if args.data_folder and args.data_folder[-1] != "/":
        args.data_folder = args.data_folder + "/"
        
    if args.results_folder and args.results_folder[-1] != "/":
        args.results_folder = args.results_folder + "/"
        
    if args.log_folder and args.log_folder[-1] != "/":
        args.log_folder = args.log_folder + "/"
    
    if args.cpus > 32:  # don't go over the amount available on the cluster
        args.cpus = 32

    args.assembler = "ray"   # TODO add more assemblers and remove this statement
    
    header = "#! /bin/sh \n"
    header += "# -*- coding: utf-8 -*- \n"
    header += "#SBATCH --cpus-per-task "+str(args.cpus)+"\n"
    header += "#SBATCH -p highmem\n"
    header += "#SBATCH --mem 245G\n"  # this is the max that worked for me
    
    plot = ""
    if args.plot:
        plot = "p"
    
    read_re = re.compile("(.*\/)*(\w*)\..*")
    
    # do the requested grouping and assembly
    assembly_names = list()
    if args.group_ids:  # custom grouping
        if len(args.group_ids) != len(args.samples):
            raise BaseException("number of group ids does not match the number of input files")
        regestry = dict()
        for i in range(len(args.group_ids)):
            try:
                regestry[args.group_ids[i]]
            except KeyError:
                regestry[args.group_ids[i]] = list()
            regestry[args.group_ids[i]].append(args.samples[i])  # dictionary of group labels with list of samples belonging to it
        for group in regestry.keys():
            assembly_names.append(group)
            assemble(regestry[group], group, args.script_folder, args.data_folder, args.results_folder, args.log_folder, header, assembler=args.assembler, cpus=args.cpus)
    elif args.pool:  # group all together
        name = "pooled"
        assembly_names.append(name)
        r = list()  # input must be a list
        for s in args.samples:  # in case of single input
            r.append(s)
        assemble(r, name, args.script_folder, args.data_folder, args.results_folder, args.log_folder, header, assembler=args.assembler, cpus=args.cpus)
    elif args.separate:  # no grouping
        #TODO don't overwrite existing .sh files
        for read in args.samples:
            read_name = read_re.match(read)
            if read_name:
                r = list()  # input must be a list
                r.append(read)
                name = read_name.group(2)
                assembly_names.append(name)
                assemble(r, name, args.script_folder, args.data_folder, args.results_folder, args.log_folder, header, assembler=args.assembler, cpus=args.cpus)
            else:
                print("invalid file name: " + str(read))
    else:  # automatic grouping
        # get unique sample names
        names = set()
        for read in args.samples:
            read_name = read_re.match(read)
            if read_name:
                name = read_name.group(2)
                names.add(name)
        # pull out and assemble all the files associated with each sample name
        for n in names:
            group = list()
            for f in args.samples:
                if f.find(n) >= 0:
                    group.append(f)
            assembly_names.append(n)
            assemble(group, n, args.script_folder, args.data_folder, args.results_folder, args.log_folder, header, assembler=args.assembler, cpus=args.cpus)

    # analysis of all the assemblies
    ill = "parse_ale.py -o "+args.results_folder+"ale_statistics.csv"+" -s"+plot+" "
    for an in assembly_names:
        ill += args.log_folder+an+"_ale.log "

    tie = "parse_bt2.py -o "+args.results_folder+"bt2.csv "
    for an in assembly_names:
        tie += args.log_folder+an+"_bt2.log "
    
    sour = "sourmash compare -k 21 --csv "+args.results_folder+"sourmash_compare.csv "
    for an in assembly_names:
        sour += args.data_folder+an+"/"+an+".sig "
    sour += "2>"+args.log_folder+"bt2_build.log "

    sour += "\n\n"
    sour += "parse_sourmash.py -"+plot+"c "+args.results_folder+"sourmash_compare.csv -o "+args.results_folder+"sourmash_parsed 2>"+args.log_folder+"parse_sourmash.log "
    
    anal_cpus = len(assembly_names)
    if anal_cpus > 32:
        anal_cpus = 32
    meta = "metaquast.py --threads "+str(anal_cpus)+" --output-dir "+args.results_folder
    meta += " --gene-finding --no-check --no-plots --no-icarus --no-snps --max-ref-number 0 "
    meta += " 2>"+args.log_folder+"metaquast.log "
    meta += "--labels \""
    for an in assembly_names:
        meta += an+","
    meta = meta[:-1]  # removeing last comma
    meta += "\" "
    for an in assembly_names:
        meta += args.data_folder+an+"/Contigs.fasta "
    
    combine = "combine_stats.py -o "+args.results_folder+"final_stats.csv"
    combine += " -a "+args.results_folder+"ale_statistics.csv"    
    combine += " -b "+args.results_folder+"bt2.csv"
    combine += " -m "+args.results_folder+"transposed_report.tex"
    if len(assembly_names) > 1:
        combine += " --other "+args.results_folder+"sourmash_parsed.csv"
    combine += " 2>"+args.log_folder+"combine_stats.log "
    
    header = "#! /bin/sh \n"
    header += "# -*- coding: utf-8 -*- \n"
    header += "#SBATCH --cpus-per-task "+str(anal_cpus)+"\n"
    header += "#SBATCH -p batch\n"
    header += "#SBATCH --mem 10G\n\n"
    
    anal = header+ tie + "\n\n" + ill + "\n\n" 
    if len(assembly_names) > 1:
        anal += sour + "\n\n" 
    anal += meta + "\n\n" + combine
    
    with open(args.script_folder+"assembly_comparison.sh", "w") as file:
        file.write(anal)
    

def assemble(sample, assembly_name, output, assembly_folder, analysis_folder, log_folder, header, cpus, assembler="ray", k=41, min_length=200):
    """Creates and saves the shell script to assemble and analyze all the input reads"""
    assembly_name = str(assembly_name)
    forward, reverse, inter, singles = extract_files(sample)
    if not output:
        output = ""
        
    assembly_folder = assembly_folder+assembly_name+"/"
    analysis_folder = analysis_folder+assembly_name+"/"
        
    command = header + "\n\n"
    if assembler.lower() == "ray":
        command += ray_meta_prep(forward, reverse, inter, singles, assembly_folder, k, min_length)
    # TODO add other assemblers here
    command += "\n\n"
    
    command += "bowtie2-build --threads "+str(cpus)+" "+assembly_folder+"Contigs.fasta 2>"+log_folder+assembly_name+"_bt2_build.log "
    command += analysis_folder+assembly_name
    command += "\n\n"
    
    bow = ""
    f = comma_separate(forward)
    r = comma_separate(reverse)
    i = comma_separate(inter)
    s = comma_separate(singles)
    
    if f:
        bow += " -1 " + f
    if r:
        bow += " -2 " + r
    if s:
        bow += " -U " + s
    if i:
        bow += " --interleaved " + i

    command += "bowtie2 -q -p "+str(cpus)+" --very-sensitive -x "+analysis_folder+assembly_name+bow+" 2>"+log_folder+assembly_name+"_bt2.log"
    command += " | samtools sort -@ "+str(cpus)+" -l 9 -O bam -T "+assembly_name+" -o "+analysis_folder+assembly_name+".bam -"
    command += "\n\n"
    
    command += "ALE --qOff 33 --metagenome "+analysis_folder+assembly_name+".bam "
    command += assembly_folder+"Contigs.fasta "
    command += analysis_folder+assembly_name+".ale "
    command += ">/dev/null 2>"+log_folder+assembly_name+"_ale.log"
    command += "\n\n"
    command += "sourmash compute -o "+assembly_folder+assembly_name+".sig "+assembly_folder+"Contigs.fasta 2>"+log_folder+assembly_name+"_sourmash_signature.log"
    
    with open(output+assembly_name+".sh", "w") as file:
        file.write(command)

    
def comma_separate(in_list):
    """separates the contents of a list with a comma and a space"""
    s = ""
    for li in in_list:
        s += li + ", "
    if s:
        s = s[:-2]
    return s


def ray_meta_prep(forward=None, reverse=None, inter=None, singles=None, output=None, k=41, min_length=200):
    """Creates the command to run Ray meta on the input files"""
    arguments = " "
    
    if output:
        arguments += "-o "+ str(output) + " "
        
    if inter:
        arguments += "-i "
        for i in inter:
            arguments += i + " "
            
    if forward and reverse:
        if len(forward) == len(reverse):
            for i in range(len(forward)):
                arguments += "-p " + forward[i] + " " + reverse[i] + " "
        else:
            raise BaseException("the number of forward reads does not match the number of reverse reads")
            
    if singles:
        arguments += "-s "
        for i in singles:
            arguments += i + " "
            
    return "mpirun Ray Meta -k "+str(k)+" -minimum-contig-length "+str(min_length)+arguments

    
def extract_files(files):
    """Separates the different input file types"""
    assert isinstance(files, list)
    re_s = re.compile(".*\.singles\..*")
    re_i = re.compile(".*\.interleaved\..*")
    re_f = re.compile(".*\.forward\..*")
    re_r = re.compile(".*\.reverse\..*")
    
    singles = list()
    inter = list()
    forward = list()
    reverse = list()
    
    for file in files:
        is_s = re_s.match(file)
        is_i = re_i.match(file)
        is_f = re_f.match(file)
        is_r = re_r.match(file)
        
        if is_s:
            singles.append(is_s.group(0))
        if is_i:
            inter.append(is_i.group(0))
        if is_f:
            forward.append(is_f.group(0))
        if is_r:
            reverse.append(is_r.group(0))
    return forward, reverse, inter, singles


if __name__ == "__main__":
    main()
