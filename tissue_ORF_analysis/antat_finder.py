"""
From Beaver et al. find the AnTat1.1 like ORFs
Smith 2024
Author: Jaclyn Smith
"""
import os
antat_fasta = "antat.fa"

os.system("bowtie-build " + antat_fasta + " " + "antat")

with open("all_orfs.fa", "r") as input, open("all_tiles.fa", "w") as output:
    n = 20
    for line in input:
        if line.startswith(">"):
            name = line
            i = 0
        else:
            tiles = [line[i:i + n] for i in range(0, len(line), n)]
            for entry in tiles:
                name_pos = name.strip() +"."+ str(i)
                entry = entry.strip()
                if len(entry) < 20:
                    pass
                else:
                    output.write(">%s\n%s\n" % (name_pos, entry))
                    i += 1

os.system("bowtie -v 2 antat -f all_tiles.fa --al aligned_tiles.fa 2> output_tile_align.txt")

with open("aligned_tiles.fa", "r") as input:
    vsg_names = set()
    i = 0
    prev_name = ""
    vsg_name = ""
    for entry in input:
        if entry.startswith(">"):
            prev_name = vsg_name
            vsg_name = entry.strip().split(".")[0]
            if prev_name == vsg_name:
                i += 1
            # elif i >= 4:
            elif i >= 15:
                vsg_names.add(vsg_name)
                i = 0
with open("all_orfs.fa", "r") as orfs, open("possible_mosaics", "w") as output:
    first = True
    seq = False
    for item in orfs:
        if item.startswith(">"):
            if not first and seq:
                output.write(">%s\n%s\n" % (name, sequence))
                seq = False
            name = item.strip()
            if name in vsg_names:
                first = False
                seq = True
                sequence = ""
        elif seq:
            sequence = sequence + item.strip()
    if not first and seq:
        output.write(">%s\n%s\n" % (name, sequence))
    # seq = False

print("test")
