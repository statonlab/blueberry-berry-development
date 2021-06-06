
#Created on June 5th, 2021
#Author: Jiali
#Usage: According to the terminal N table to adjust gene location.
#python remove_terminalN_gff.py <N_tatble.txt> <gene annotation.gff3> <output.gff3>


import sys

file_table = sys.argv[1]
in_gff = sys.argv[2]
out_gff = sys.argv[3]

def forward_determine(content):
    if int(content[4]) > int(content[3]):
        return True
    else:
        return False

with open(in_gff) as in_file, open(file_table) as N_Table, open(out_gff,"w") as out:
    for line in in_file:
        if line.startswith("#"):
            continue
        else:
            content = line.split("\t")
            chr = content[0]
            for table_line in N_Table:
                table_content = table_line.strip("\n").split("\t")
                is_forward = forward_determine(content)
                if chr == table_content[0]:
                    if is_forward:
                        new_start = int(content[3]) - int(table_content[1])
                        new_end = int(content[4]) - int(table_content[1])
                    else:
                        new_start = int(content[3]) - int(table_content[2])
                        new_end = int(content[4]) - int(table_content[2])
                else:
                    continue
            new_gff = content[:3]+[str(new_start)]+[str(new_end)]+content[5:]
            out.write("\t".join(new_gff))
            N_Table.seek(0)
            