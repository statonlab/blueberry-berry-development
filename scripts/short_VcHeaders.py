"""
Created on Jan 18, 2021
Author: Jiali
Usage: python short_VcHeaders.py <input.fasta> <output.fasta>
This script is to edit the long mRNA names into truncated format to match the protein truncated names
"""

import sys
file_in = sys.argv[1]
file_out = sys.argv[2]

with open(file_in) as f, open(file_out,"w") as out:
    for line in f:
        if line.startswith(">"):
            header = line.split("||")
            ID = header[4].replace(".mRNA1","-mRNA-1")
            fix_ID = ID.replace("maker-","").replace("snap_masked-","").replace("augustus_masked-","")
            out.write(">"+fix_ID+"\n")
        else:
            out.write(line)