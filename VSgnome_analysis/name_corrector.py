def fix_name(input_file, output_file):
    with open(input_file, "r") as input,\
            open(output_file, "w") as output:
        for line in input:
            if line.startswith(">"):
                new_name = "".join(line.split(" "))
                new_name = "".join(new_name.split("\t"))
                output.write(new_name)
            else:
                output.write(line)


def main():
    """Execute the functions"""
    fix_name("./unknown_protein_VSGBlast/all_expressed_genomic_VSGs.fa",
             "./unknown_protein_VSGBlast/all_expressed_genomic_VSGs_nameCorrected.fa")
    fix_name("./unknown_protein_VSGBlast/TriTrypDB-66_TbruceiLister427_2018_unknown.fa",
             "./unknown_protein_VSGBlast/TriTrypDB-66_TbruceiLister427_2018_unknown_nameCorrected.fa")
    print("finished script")
    fix_name("TriTrypDB-66_TbruceiLister427_2018_VSGs.fa",
             "TriTrypDB-66_TbruceiLister427_2018_VSGs_nameCorrected.fa")


if __name__ == '__main__':
    main()
