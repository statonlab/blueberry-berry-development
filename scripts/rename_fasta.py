"""
created on May 22, 2020
Author: Jiali
This scirpt will change sequence headers based on the file provided
"""
import argparse

parser=argparse.ArgumentParser(description="program that replaces fasta records")
parser.add_argument("-i", help="input fasta", type=argparse.FileType('r'))
parser.add_argument("-r", help="replacement records file", type=argparse.FileType('r'))
parser.add_argument("-o", help="output file")
args = parser.parse_args()
newfasta=open(args.o,'w') 

for line in args.i:
    if line.startswith('>'):
        name = line.rstrip().strip(">")
        for matches in args.r:
            if name in matches:
                replacement = matches.strip("\n").split("\t")
                newname = replacement[1]
        newfasta.write(">"+newname+"\n")
        args.r.seek(0)
    else: 
        newfasta.write(line)