def gff_sequence_puller(VSG_gff_input_file, associated_genome_file, vsg_fasta_output):
    with open(VSG_gff_input_file, "r") as gff,\
            open(associated_genome_file, "r") as genome, \
            open(vsg_fasta_output, "w") as output:

        genome_dict = {}
        start = True
        set_desc = set()
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

        for line in gff:
            if line.startswith("##"):
                pass
            else:
                # obtain sequences based off gff coordinates. add 300bps up and downstream to fasta
                vsg_name = ">" + line
                chr = line.split("\t")[0]
                start = int(line.split("\t")[3]) - 300
                end = int(line.split("\t")[4]) + 300
                if start < 0:
                    start = 0
                elif end < 0:
                    end = 0
                vsg = genome_dict[chr][start:end]
                output.write(vsg_name)
                output.write(vsg + "\n")
                description = line.split("description=")[1].split(";")[0]
                set_desc.add(description)


def main():
    """Execute the functions"""
    gff_sequence_puller("TriTrypDB-66_TbruceiLister427_2018_VSGs.gff",
                        "./genome_2018_associatedFiles/TriTrypDB-66_TbruceiLister427_2018_Genome.fasta",
                        "TriTrypDB-66_TbruceiLister427_2018_VSGs.fa")
    gff_sequence_puller("./genomes/TriTrypDB-66_TbruceiEATRO1125_VSGs.gff",
                        "./genomes/TriTrypDB-66_TbruceiEATRO1125_Genome.fasta",
                        "./genomes/TriTrypDB-66_TbruceiEATRO1125_VSGs.fa")
    gff_sequence_puller("./genomes/TriTrypDB-66_TbruceiLister427_VSGs.gff",
                        "./genomes/TriTrypDB-66_TbruceiLister427_Genome.fasta",
                        "./genomes/TriTrypDB-66_TbruceiLister427_VSGs.fa")
    gff_sequence_puller("./genomes/TriTrypDB-66_TbruceiTREU927_VSGs.gff",
                        "./genomes/TriTrypDB-66_TbruceiTREU927_Genome.fasta",
                        "./genomes/TriTrypDB-66_TbruceiTREU927_VSGs.fa")
    gff_sequence_puller("TriTrypDB-66_TbruceiLister427_2018_unknown.gff",
                        "./genome_2018_associatedFiles/TriTrypDB-66_TbruceiLister427_2018_Genome.fasta",
                        "TriTrypDB-66_TbruceiLister427_2018_unknown.fa")

    print("finished script")


if __name__ == '__main__':
    main()
