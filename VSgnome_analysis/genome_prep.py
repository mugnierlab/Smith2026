def genome_coords(input_file, output_tsv):
    """
    inputs a genome FASTA file and outputs coordinates for all entries
    :param input_file: string, path to input file usually a genome FASTA file
    :param output_tsv: string, path to output tsv to be printed with
                               chromosome name and DNA length
    :return: output written to file, no return
    """
    with open(input_file, "r") as genome, \
            open(output_tsv, "w") as output:
        genome_dict = {}
        start = True
        for line in genome:
            if line.startswith(">"):
                if not start:
                    genome_dict[chr_name] = sequence
                chr_name = line.strip().split(" ")[0].split(">")[1]
                sequence = ""
                start = False
            else:
                sequence = sequence + line.strip()
        genome_dict[chr_name] = sequence

        for chromosome in genome_dict :
            output.write(chromosome + "\t" + str(len(genome_dict[chromosome])) + "\n")


def gff_filter(gff_unfiltered, tsv_subset, gff_output):
    with open(gff_unfiltered, "r") as gff, \
            open(tsv_subset, "r") as blast_filter, \
            open(gff_output, "w") as gff_out:
        gene_list = []
        for line in blast_filter:
            if line.startswith("unknown"):
                pass
            else:
                gene_list.append(line.split("ID=")[1].split(";")[0])
        for entry in gff:
            if entry.startswith("#"):
                gff_out.write(entry)
            else:
                gene_entry = entry.split("ID=")[1].split(";")[0]
                if gene_entry in gene_list:
                    gff_out.write(entry)


def add_duplicates_back(filtered_fasta_file, fasta_2018, duplicate_fasta):
    with open(filtered_fasta_file, "r") as input_fasta,\
         open(fasta_2018, "r") as checkList_fasta, \
         open(duplicate_fasta, "w") as duplicates:

        vsg_id_list = []
        write_next = False
        for line in input_fasta:
            if "ID" in line:
                vsg_id = line.split("ID=")[1].split(";")[0]
                vsg_id_list.append(vsg_id)

        for entry in checkList_fasta:
            if "ID" in entry:
                vsg_id = entry.split("ID=")[1].split(";")[0]
                if vsg_id not in vsg_id_list:
                    duplicates.write(entry)
                    write_next = True
            elif write_next:
                duplicates.write(entry)
                write_next = False


def check_2018_selected_as_master(clstr_file):
    with open(clstr_file, "r") as clstr:
        leader = ""
        first = True
        for line in clstr:

            if line.startswith(">") and not first:
                if "ID=" in leader:
                    pass
                    # print("leader is a 2018 VSG")
                elif "ID=" not in leader:
                    for entry in cluster:
                        if "ID=" in entry:
                            # print("leader is not 2018 VSG and it is a group member")
                            break
                        else:
                            pass
                        # print("group is just from VSGnome")
                    #print()
                cluster = []

            elif not first:
                cluster.append(line.split(", ")[1])
                if "*" in line:
                    leader = line.split(", ")[1]
            else:
                first = False
                cluster = []



def main():
    """Execute the functions"""


    # genome_coords("./genome_2018_associatedFiles/TriTrypDB-66_TbruceiLister427_2018_Genome.fasta",
    #               "genome_2018_coords.tsv")
    # gff_filter("./TriTrypDB-66_TbruceiLister427_2018_unknown.gff",
    #            "./unknown_protein_VSGBlast/blastResults_filtered.tsv",
    #            "./TriTrypDB-66_TbruceiLister427_2018_unknown_filtered.gff")
    #
    # add_duplicates_back("all_unique_Lister427VSGs.fa",
    #                     "TriTrypDB-66_TbruceiLister427_2018_VSGs_nameCorrected.fa",
    #                     "TriTrypDB-66_TbruceiLister427_2018_VSGs_nameCorrected_duplicates.fa")
    #check_2018_selected_as_master("all_unique_Lister427VSGs.fa.clstr")
    add_duplicates_back("all_unique_Lister427VSGs.fa",
                                    "TriTrypDB-66_TbruceiLister427_2018_unknown_filtered_nameCorrected.fa",
                                    "TriTrypDB-66_TbruceiLister427_2018_unknown_filtered_nameCorrected_duplicates.fa")

    print("finished script")


if __name__ == '__main__':
    main()
