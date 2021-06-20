"""
Created on June 16th, 2021
Author: Jiali
Usage: gFACS removing the CDS information. According to the gFACS filtering results, keep the genes that retained in gFACS from the original braker output.
python filter_braker.py <input protein fasta> <input braker.gff3> <output gff3>
"""

import sys

fasta_file = sys.argv[1]
file_in = sys.argv[2]
file_out = sys.argv[3]

gene_list = [] # make an empty list, put the protein IDs in the list
with open(fasta_file) as f:
    for line1 in f:
        if line1.startswith(">"):
            ID = line1.replace(">","")
            gene_list.append(ID)

with open(file_in) as f_in, open(file_out, "w") as out:
    for line2 in f_in:
        content = line2.split("\t")
        info = content[8].split(".")
        if info[0] in gene_list:
            out.write(line2)