"""
Created on Jan 18th, 2021
Author: Jiali
Usage: python grepGeneIDs.py <input.header>
"""

import sys
file_in = sys.argv[1]

with open(file_in) as f:
    for line in f:
        content=line.split("||")
        ID = content[4].replace(".mRNA1","-mRNA-1")
        fix_ID = ID.replace("maker-","").replace("snap_masked-","").replace("augustus_masked-","")
        print(fix_ID)