"""
Create on March 5, 2012
Author: Jiali
Usage: edit_hap2_chrs.py <input.fasta> <output.fasta>
"""
import sys
file_in = sys.argv[1]
file_out = sys.argv[2]

with open(file_in, "r") as f, open(file_out,"w") as out:
    for line in f:
        if line.startswith(">chr"):
            header = line.strip(">")
            chr_num = int(header.replace("chr",""))
            new_chr_num = chr_num - 12
            new_chr = "Chr"+str(new_chr_num)+"_H2"
            out.write(">"+new_chr+"\n")
        else:
            out.write(line)
