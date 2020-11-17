"""
This program is to get the GO IDs for blueberry genes, from the entap output file. 
Usage: python extract_GO.py <entap.table> <output.table>
"""

import sys
file_in = sys.argv[1]
file_out = sys.argv[2]

with open(file_in, "r") as f, open(file_out,"w") as out:
    for line in f:
        if "GO" in line:
            content = line.split("\t")
            out.write(content[0]+ "\t" + '\t'.join(content[-7:-2])+"\n")