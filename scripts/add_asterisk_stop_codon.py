"""
Usage python add_asterisk_stop_codon.py <input fasta> <output fasta>
"""
import sys
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

fasta_in = sys.argv[1]
fasta_out = sys.argv[2]

def add_stop_codon(record):
    return SeqRecord(seq=record.seq+"*", \
            id=record.id, \
            description="")

fasta_seq = SeqIO.parse(open(fasta_in), "fasta")
with open(fasta_out, "w") as out:
    for seq in fasta_seq:
        seq_with_stop = add_stop_codon(seq)
        SeqIO.write(seq_with_stop, out, "fasta")