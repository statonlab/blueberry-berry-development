"""
This program take the results from synmap to generate the correct format for circos
"""

import sys
file_in = sys.argv[1]
output = sys.argv[2]

with open(file_in) as f, open(output,"w") as out:
    for line in f:
        if line.startswith("#"):
            continue
        else:
            genome2 = line.strip("\n").split("\t")[1]
            genome1 = line.strip("\n").split("\t")[5]
            first_coor = genome1.split("||")
            second_coor = genome2.split("||")
            identity = float(first_coor[8])
            if identity > 90.0:
                #color = first_coor[0].replace("HiC_scaffold_", "color=chr")
                color = "color=chr"+second_coor[0]
                out.write("chr"+first_coor[0] + "\t" + first_coor[1] + "\t" + first_coor[2] + "\t" +
                "chr"+second_coor[0] + "\t" + second_coor[1]+ "\t" + second_coor[2] + "\t"+
                color + "\n")

