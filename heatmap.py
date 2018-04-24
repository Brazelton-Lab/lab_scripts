#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Copyright:

    heatmap  Produces a heatmap from the raw sourmash output.

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
__version__ = '1.0'


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn
import argparse


def main():
    ap = argparse.ArgumentParser(description="Produces a heatmap from the raw sourmash output.")
    ap.add_argument("input", help="The name of the input file. It typically lacks a file type ending.")
    ap.add_argument("lables", help="A file containing an ordered list of lables on each line corresponding to the files compared by sourmash.")
    ap.add_argument("-s", "--suffix", default="_clustermap.png", help="The end of the file name of the output plot files.")
    args = ap.parse_args()
    
    exp = args.input
    labels = args.lables
    names = []
    d = np.load(exp)
    with open(labels) as file:
        for line in file:
            if line:
                names.append(line)
    data = pd.DataFrame(d, index=names, columns=names)
    
    cm = seaborn.clustermap(data)
    plt.setp(cm.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)		
    cm.savefig(exp + args.sufix)


if __name__ == "__main__":
    main()

