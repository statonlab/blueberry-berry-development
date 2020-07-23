"""
Created on July 21, 2020
Author: Jiali
This script is used to flip the sequences
Usage: python fix_scaff13.py <input_fasta.fa> <int> <output.fasta>
"""
import argparse

parser=argparse.ArgumentParser(description="takes a location to break the sequences and flip the two parts")
parser.add_argument("-i", help="input fasta")
parser.add_argument("-loc", help="break location [integer]", type=int)
parser.add_argument("-o", help="output file")
args = parser.parse_args()
newfasta=open(args.o,'w')

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# read the fasta file and get the sequence
fasta = SeqIO.parse(open(args.i), "fasta")
for fasta_seq in fasta:
    seq = fasta_seq.seq
    first_part = seq[0:args.loc]
    second_part = seq[args.loc:]
    new_seq = second_part + first_part

# save the header and new sequences into a seq record
    new_record = SeqRecord(seq = new_seq, id = fasta_seq.id, description="")
    SeqIO.write(new_record, newfasta, "fasta")