"""
Created on June 1st, 2021
Author: Jiali 
Usage: This program is to remove the Ns in the begining and the end of the scaffolds/contigs, and output a table how many Ns being removed.
python remove_terminalN.py <input fasta> <output fasta> <output table txt>
"""

import sys
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

file_in = sys.argv[1]
file_out = sys.argv[2]
table_out = sys.argv[3]

def make_seq_record(record, new_Seq):
    return SeqRecord(seq=new_Seq, \
            id=record.id, \
            description="")

fasta_seq = SeqIO.parse(open(file_in), "fasta")
with open(file_out, "w") as out, open(table_out,"w") as table:
    for seq in fasta_seq:
        header = seq.id
        sequence = seq.seq
        begin = sequence[:30]
        begin_N_num = begin.count("N")
        end = sequence[-30::]
        end_N_num = end.count("N")
        new_Seq = begin.lstrip("N")+sequence[30:-30]+end.rstrip("N")
        out_record = make_seq_record(seq,new_Seq) 
        SeqIO.write(out_record, out, "fasta")
        table.write(header+"\t"+str(begin_N_num)+"\t"+str(end_N_num)+"\n")