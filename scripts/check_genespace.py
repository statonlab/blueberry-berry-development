"""
Created on June, 1st 2021
Author: Jiali
Usage: Finding if the gene content in the gff is in the one direction.
python check_genespace.py <input gff file>
"""

import sys

file_in = sys.argv[1]

forward = 0
reverse = 0

with open(file_in) as f:
    for line in f:
        if line.startswith("#"):
            continue
        else:
            content = line.split("\t")
            start = content[3]
            end = content[4]
            if start < end:
                forward = forward+1
            else:
                reverse = reverse+1

    print("forward direction: "+str(forward))
    print("reverse direction: "+str(reverse))