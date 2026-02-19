"""
This script is for finding duplicate VSG families which do not have any diversification
"""

import glob
import itertools

cluster_dict = {}
for file_name in glob.glob("./cluster_BLAST/Lister*"):
    with open(file_name, "r") as input_file:
        for line in input_file:
            if line.startswith("VSG"):
                continue
            try:
                cluster_dict[line.split(",")[1]].append(line.split(",")[0])
            except KeyError:
                cluster_dict[line.split(",")[1]] = [line.split(",")[0]]

    vsg_dict = {}
    if "Only" in file_name:
        path = "./vsgnomes_2018_renamed/all_unique_posDuplicate_Lister427VSGsonly_AnTat_renamed.fa"
    else:
        path = "./vsgnomes_2018_renamed/all_unique_posDuplicate_Lister427VSGs_AnTat_renamed.fa"
    with open(path, "r") as vsg_seqs:
        notStart = False
        name = ""
        seq = ""
        for line in vsg_seqs:
            if line.startswith(">"):
                if notStart:
                    vsg_dict[name] = seq
                    seq = ""
                name = line.strip().split(">")[1].split(" : ")[0].split(" GB ")[0]
                notStart = True
            else:
                seq += line.strip()

    vsg_dict["AnTat1.1"] = seq

    duplicated_vsgs = []

    for item in cluster_dict:
        same = False
        overwrite = False
        if len(cluster_dict[item]) > 1:
            for entry in itertools.product(cluster_dict[item], cluster_dict[item]):
                if vsg_dict[entry[0]] != vsg_dict[entry[1]]:
                    print(entry)
                    same = False
                    overwrite = True
                    break
                elif entry[0] != entry[1]:
                    same = True
                    print("these are the same")
            if not overwrite and same:
                duplicated_vsgs = duplicated_vsgs + cluster_dict[item]
    print("final")
