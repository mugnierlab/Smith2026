def VSG_puller(gff_file_path, gff_file_output):
    with open(gff_file_path, "r") as gff, \
            open(gff_file_output, "w") as output:
        for line in gff:
            if line.startswith("##"):
                output.write(line)
            else:
                # remove start with chromosome specification, sometimes includes "VSG"
                line_list = line.split("\t")
                line_list.pop(0)
                test_line = "".join(line_list)
                if "VSG" in test_line and "protein_coding_gene" in test_line and "Exclusion" not in test_line \
                        and "exclusion" not in test_line:
                    output.write(line)
                elif "VSG" in test_line and "pseudogene" in test_line and "pseudogenic_transcript" not in test_line \
                        and "Exclusion" not in test_line and "exclusion" not in test_line:
                    output.write(line)
                # removed extra genes pulled by variant key word. multiple versions of invariant and RING-variant proteins
                elif "variant" in test_line and "invariant" not in test_line and "Invariant" not in test_line \
                        and "RING" not in test_line and "protein_coding_gene" in test_line and "Histone" not in test_line and "histone" not in test_line \
                        and "ubiquitin" not in test_line and "Exclusion" not in test_line and "Nonvariant" not in test_line \
                        and "nonvariant" not in test_line:
                    output.write(line)
                elif "variant" in test_line and "invariant" not in test_line and "Invariant" not in test_line \
                        and "RING" not in test_line and "pseudogene" in test_line \
                        and "pseudogenic_transcript" not in test_line and "Histone" not in test_line and "histone" not in test_line \
                        and "ubiquitin" not in test_line and "Exclusion" not in test_line and "Nonvariant" not in test_line \
                        and "nonvariant" not in test_line:
                    output.write(line)


def pull_unknown_genes(gff_file_path, gff_file_output):
    with open(gff_file_path, "r") as gff, \
            open(gff_file_output, "w") as output:
        for line in gff:
            if line.startswith("##"):
                output.write(line)
            else:
                # remove start with chromosome specification, sometimes includes "VSG"
                line_list = line.split("\t")
                line_list.pop(0)
                test_line = "".join(line_list)
                if "VSG" in test_line:
                    pass            # removed extra genes pulled by variant key word. multiple versions of invariant and RING-variant proteins
                elif "variant" in test_line and "invariant" not in test_line and "Invariant" not in test_line \
                        and "RING" not in test_line and "Histone" not in test_line and "histone" not in test_line \
                        and "ubiquitin" not in test_line and "Exclusion" not in test_line and \
                        "exclusion" not in test_line and "Nonvariant" not in test_line and \
                        "nonvariant" not in test_line:
                    pass
                elif "hypothetical protein" in test_line and "pseudogene" in test_line and "Parent" not in test_line:
                    output.write(line)
                elif "hypothetical protein" in test_line and "protein_coding_gene" in test_line and "Parent" not in test_line:
                    output.write(line)
                elif "pseudogene" in test_line and "pseudogenic_transcript" not in test_line \
                        and "Parent" not in test_line:
                    output.write(line)
                elif "pseudogene" in test_line and "protein_coding_gene" in test_line and "Parent" not in test_line:
                    output.write(line)
                elif "unknown" in test_line and "pseudogenic_transcript" not in test_line \
                        and "Parent" not in test_line:
                    output.write(line)
                elif "unknown" in test_line and "protein_coding_gene" in test_line and "Parent" not in test_line:
                    output.write(line)



def main():
    """Execute the functions"""
    VSG_puller("./genome_2018_associatedFiles/TriTrypDB-66_TbruceiLister427_2018.gff",
               "TriTrypDB-66_TbruceiLister427_2018_VSGs.gff")
    VSG_puller("./genomes/TriTrypDB-66_TbruceiEATRO1125.gff",
               "./genomes/TriTrypDB-66_TbruceiEATRO1125_VSGs.gff")
    VSG_puller("./genomes/TriTrypDB-66_TbruceiLister427.gff",
               "./genomes/TriTrypDB-66_TbruceiLister427_VSGs.gff")
    VSG_puller("./genomes/TriTrypDB-66_TbruceiTREU927.gff",
               "./genomes/TriTrypDB-66_TbruceiTREU927_VSGs.gff")
    pull_unknown_genes("./genome_2018_associatedFiles/TriTrypDB-66_TbruceiLister427_2018.gff",
                       "TriTrypDB-66_TbruceiLister427_2018_unknown.gff")
    print("finished script")


if __name__ == '__main__':
    main()
