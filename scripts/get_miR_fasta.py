"""
Created on June 2021
Author: Jiali
Usage: extract miRNA sequence from csv and save into fasta 
"""

import sys

file_in = sys.argv[1]
file_out = sys.argv[2]

with open(file_in) as f, open(file_out,"w") as out:
    line1 = f.readline()
    for line in f:
        content = line.split(",")
        ID = content[0]
        seq = content[12]
        out.write(">"+ID+"\n"+seq+"\n")
