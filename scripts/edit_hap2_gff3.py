"""
Create on March 5, 2012
Author: Jiali
Usage: edit_hap2_gff3.py <input.gff3> <output.gff3>
"""
import sys
file_in = sys.argv[1]
file_out = sys.argv[2]

with open(file_in, "r") as f, open(file_out,"w") as out:
    for line in f:
        if line.startswith("chr"):
            content = line.split("\t")
            chr_num = int(content[0].replace("chr",""))
            new_chr_num = chr_num - 12
            new_chr = "Chr"+str(new_chr_num)+"_H2"
            new_content = [new_chr]+content[1:]
            out.write("\t".join(new_content))
        else:
            out.write(line)