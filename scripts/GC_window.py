"""
Created on Dec 9th
Author: Jiali
Calculate GC content within a sliding window difined by user

Usage:
python GC_window.py --input <input fasta> --output <output file name> --window <int>
"""
import argparse
import sys
from Bio import SeqIO

parser = argparse.ArgumentParser(description="CG content window analysis", usage="%(prog)s --input <input file> --output <output file> --window <int>")
parser.add_argument("--input", type=str, help="input fasta file")
parser.add_argument("--output", type=str, help="output bedgraph file")
parser.add_argument('--window', action='store', type=int, required=True, help="sliding window size")
args = parser.parse_args()
fafile = args.input
output_file = args.output
window = args.window

# calculate GC content of a given sequence in percentage
def computeGC(seq):
    sequence = seq.upper()
    C_num = sequence.count("C")
    G_num = sequence.count("G")
    GC = (C_num + G_num)/len(sequence)
    return round(GC * 100, 1)

# read fasta file
fasta_seq = SeqIO.parse(open(fafile), "fasta")

# create a list of number to slide fasta sequence
def sliding(window, chro_size):
    if window < chro_size:
        starts = int(chro_size/window)
        return [x * window for x in range(0, starts)]
    else:
        return None

if __name__ == "__main__":
    with open(output_file, "w") as out:
        for seq in fasta_seq:
            size = len(seq)
            starts = sliding(window, size)
            if starts == None:
                continue
            for i, start in enumerate(starts):
                if i < len(starts) - 1:
                    segment = seq.seq[start: start + window]
                    GC_content = computeGC(segment)
                    out.write(seq.id + "\t"+str(start)+"\t"+str(start + window)+"\t"+str(GC_content)+"\n")
                else:
                    segment = seq.seq[start: size-1]
                    GC_content = computeGC(segment)
                    out.write(seq.id + "\t"+str(start)+"\t"+str(size)+"\t"+str(GC_content)+"\n")

