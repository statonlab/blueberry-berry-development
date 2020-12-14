"""
Created on Dec 14, 2020
Author: Jiali
This program edit the density table from DensityMap.pl and create the format as bedgraph for circos to plot heatmap
Usage: edit_density.py <input.csv> <output.bg>
"""

import sys
file_in = sys.argv[1]
file_out = sys.argv[2]

with open(file_in) as f, open(file_out,"w") as out:
    next(f) # skip the header line
    for line in f:
        content = line.split("\t")
        content.pop(1) # remove the second element which is the feature name
        out.write("\t".join(content))