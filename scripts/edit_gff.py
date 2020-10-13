"""
Created on Oct 12th, 2020
Author: Jiali
This script is to format the gff3 files into the way DupGen_finder input wanted.
Usage: python edit_gff.py <input_gff> <output_gff>
"""
import sys

file_in = sys.argv[1]
file_out = sys.argv[2]

with open(file_in) as f, open(file_out,"w") as out:
    for line in f:
        line_content = line.strip("\n").split("\t")
        gene_names = line_content[8].split(";")
        gene_id = gene_names[1].replace("Name=","")
        reformat = [line_content[0]] + [gene_id] + line_content[3:5]
        reformat_str = "\t".join(reformat)
        out.write(reformat_str+"\n")
