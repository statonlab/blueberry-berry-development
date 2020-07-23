"""
Usage python reverse_complement.py <input fasta> <input headers> <output fasta>
"""
import sys
from Bio.Seq import Seq
from Bio import SeqIO

fasta_in = sys.argv[1]
tab_file = sys.argv[2]
fasta_out = sys.argv[3]

contigs_list = []
with open(tab_file, "r") as contigs:
    for line in contigs:
        header = line.strip("\n")
        contigs_list.append(header)

wanted = set(contigs_list)
print(contigs_list)

def make_rc_record(record):
    return SeqRecord(seq=record.seq.reverse_complement(), \
            id=record.id, \
            description="reverse complement")

fasta_seq = SeqIO.parse(open(fasta_in), "fasta")
with open(fasta_out, "w") as out:
    for seq in fasta_seq:
        if seq.id in wanted:
            seq_rc = make_rc_record(seq)
            print(seq.id)
            SeqIO.write(seq_rc, out, "fasta")
        else:
            SeqIO.write([seq], out, "fasta")