"""
change the format of gff3 to extract the exon sequences
"""

import sys

file_in = sys.argv[1]
file_out = sys.argv[2]

with open(file_in) as f, open(file_out, "w") as out:
    for line in f:
        content = line.split("\t")
        ID = content[8].split(";")
        content[2] = ID[0].replace("ID=","")
        out.write("\t".join(content))
