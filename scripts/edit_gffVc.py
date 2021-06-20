"""
Created on Jan 18th, 2021
Author: Jiali
This script is to format the Vc bed files into the way DupGen_finder input wanted.
Usage: python edit_gff.py <input_gff> <output_gff>
"""
import sys

file_in = sys.argv[1]
file_out = sys.argv[2]

with open(file_in) as f, open(file_out,"w") as out:
    for line in f:
        if "gene" in line:
            line_content = line.strip("\n").split("\t")
            gene_names = line_content[3]+"-mRNA-1"
            fix_ID = gene_names.replace("maker-","").replace("snap_masked-","").replace("augustus_masked-","")
            reformat = [line_content[0]] + [fix_ID] + line_content[1:3]
            reformat_str = "\t".join(reformat)
            out.write(reformat_str+"\n")