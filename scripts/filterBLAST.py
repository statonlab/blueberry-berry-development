"""
Created on Aug 21st, 2019
Author: Jiali

This program takes the tab seperated BLASTn results as input, calculating the cumulative identity percentage (CIP)
and cumulative alignment length percentage (CALP) as two parameters to filtering out the blast hits less than 
60% CIP and less than 70% CALP. It produces the highest cumulative percentage identity over the longest cumulative 
length thereby increasing stringency in defining conservation between two genome sequences. (Jerome Salse, 2009)

Usage:
python filterBLAST.py <input file name> <output file name> --cip int --calp int
"""
import argparse
import sys
parser = argparse.ArgumentParser(description="Filter blast hits based on CIP and CALP", usage="%(prog)s <input file> <output file> [--cip int] [--calp int]")
parser.add_argument("input", type=str, help="input blast tabular file")
parser.add_argument("output", type=str, help="output file")
parser.add_argument('--cip', action='store', type=int, required=True, help="cumulative identity percentage cutoff")
parser.add_argument('--calp', action='store', type=int, required=True, help="cumulative alignment length percentage cutoff")
args = parser.parse_args()
BLASTfile = args.input
output_file = args.output

def is_overlap(a, b):
    if not (int(a[0]) > int(b[1]) or int(a[1]) < int(b[0])):
        return True

def LongestPara(HSP_length, HSP_iden):
    max_length = max(HSP_length)
    max_iden = HSP_iden[HSP_length.index(max_length)]
    return [max_length], [max_iden]        

def getHSP(ID):
    HSP_length=[]
    HSP_iden=[]
    hits = []
    with open(BLASTfile) as infile:
        for line in infile:
            if ID in line:
                hit_info = line.strip("\n").split("\t")
                gaps = int(hit_info[6])
                align_len = int(hit_info[4]) - gaps
                HSP_length.append(align_len)
                HSP_iden.append(hit_info[3])
                hits.append(hit_info[7:9])
                query_len = hit_info[2]
        return HSP_iden, HSP_length, query_len, hits

def removeOverlap(HSP_iden, HSP_length, hits):
    non_overlap_length = []
    non_overlap_iden =[]
    overlap_length = []
    overlap_iden = []
    a_hit = hits[0]
    for i, hit in enumerate(hits):
        overlapped = is_overlap(a_hit, hit)
        if overlapped:
            overlap_length.append(HSP_length[i])
            overlap_iden.append(HSP_iden[i])
        else:
            non_overlap_length.append(HSP_length[i])
            non_overlap_iden.append(HSP_iden[i])
            a_hit = hit
    max_overlap_length, max_overlap_iden = LongestPara(overlap_length, overlap_iden)
    all_length = max_overlap_length + non_overlap_length
    all_iden = max_overlap_iden + non_overlap_iden
    return all_iden, all_length

def getParameters(HSP_iden, HSP_length, query_len):
    AL = sum(HSP_length)
    CIP = 0.0
    for i, iden in enumerate(HSP_iden):
        IP = float(iden)*(float(HSP_length[i])/float(AL))
        CIP = CIP+IP
    CALP = (float(AL)/float(query_len)) * 100
    return CIP, CALP

def main():
    with open(BLASTfile) as file_in, open(output_file, "w") as file_out:
        for line in file_in:
            ID_lines = line.strip("\n").split("\t")
            HSPair = "\t".join(ID_lines[0:2])
            HSP_iden, HSP_length, query_len, hits = getHSP(HSPair)
            non_overlap_iden, non_overlap_length  = removeOverlap(HSP_iden, HSP_length, hits)
            CIP, CALP = getParameters(non_overlap_iden, non_overlap_length, query_len)
            #print round(CIP,2), round(CALP,2)
            if CIP > args.cip and CALP > args.calp:
                file_out.write(line)

main()