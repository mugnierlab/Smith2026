"""This script runs all necessary functions to generate a more comprehensive VSG
   annotation for the 2018 genome"""

import vsg_sequence_pull
import genome_prep
import gff_vsg_pull
import name_corrector
import os

blast_folder = "unknown_VSG_Blast"


def obtain_blast_tsv():
    try:
        os.mkdir(blast_folder)
        print("folder '{}' created ".format(blast_folder))
    except FileExistsError:
        print("folder {} already exists".format(blast_folder))
    # first filter gffs for VSG containing genes. these have been checked to ensure only VSGs
    # were pulled out via confirming the descriptions by eye
    # this also includes flanking sequences +/- 300 bps in each end.
    gff_vsg_pull.VSG_puller("./genome_2018_associatedFiles/TriTrypDB-66_TbruceiLister427_2018.gff",
                            "TriTrypDB-66_TbruceiLister427_2018_VSGs.gff")
    gff_vsg_pull.VSG_puller("./genomes/TriTrypDB-66_TbruceiEATRO1125.gff",
                            "./genomes/TriTrypDB-66_TbruceiEATRO1125_VSGs.gff")
    gff_vsg_pull.VSG_puller("./genomes/TriTrypDB-66_TbruceiLister427.gff",
                            "./genomes/TriTrypDB-66_TbruceiLister427_VSGs.gff")
    gff_vsg_pull.VSG_puller("./genomes/TriTrypDB-66_TbruceiTREU927.gff",
                            "./genomes/TriTrypDB-66_TbruceiTREU927_VSGs.gff")

    # find all the unknown genes in the 2018 genome with flanking sequences
    gff_vsg_pull.pull_unknown_genes("./genome_2018_associatedFiles/TriTrypDB-66_TbruceiLister427_2018.gff",
                                    "TriTrypDB-66_TbruceiLister427_2018_unknown.gff")

    # get a FASTA associated with the gff file
    # works for any gff file, genome and output file
    vsg_sequence_pull.gff_sequence_puller("TriTrypDB-66_TbruceiLister427_2018_VSGs.gff",
                                          "./genome_2018_associatedFiles/" +
                                          "TriTrypDB-66_TbruceiLister427_2018_Genome.fasta",
                                          "TriTrypDB-66_TbruceiLister427_2018_VSGs.fa")
    vsg_sequence_pull.gff_sequence_puller("./genomes/TriTrypDB-66_TbruceiEATRO1125_VSGs.gff",
                                          "./genomes/TriTrypDB-66_TbruceiEATRO1125_Genome.fasta",
                                          "./genomes/TriTrypDB-66_TbruceiEATRO1125_VSGs.fa")
    vsg_sequence_pull.gff_sequence_puller("./genomes/TriTrypDB-66_TbruceiLister427_VSGs.gff",
                                          "./genomes/TriTrypDB-66_TbruceiLister427_Genome.fasta",
                                          "./genomes/TriTrypDB-66_TbruceiLister427_VSGs.fa")
    vsg_sequence_pull.gff_sequence_puller("./genomes/TriTrypDB-66_TbruceiTREU927_VSGs.gff",
                                          "./genomes/TriTrypDB-66_TbruceiTREU927_Genome.fasta",
                                          "./genomes/TriTrypDB-66_TbruceiTREU927_VSGs.fa")
    vsg_sequence_pull.gff_sequence_puller("TriTrypDB-66_TbruceiLister427_2018_unknown.gff",
                                          "./genome_2018_associatedFiles/" +
                                          "TriTrypDB-66_TbruceiLister427_2018_Genome.fasta",
                                          "./" + blast_folder +
                                          "/TriTrypDB-66_TbruceiLister427_2018_unknown.fa")

    os.system("cat ./genomes/*VSGs.fa abex7_all_orfs.fa > ./genomes/all_genomic_annotated_VSGs.fa")
    os.system("cat ./genomes/all_genomic_annotated_VSGs.fa ./cross_VSGnomes/* >" +
              " ./" + blast_folder + "/all_expressed_genomic_VSGs.fa")

    name_corrector.fix_name("./" + blast_folder + "/all_expressed_genomic_VSGs.fa",
                            "./" + blast_folder + "/all_expressed_genomic_VSGs_nameCorrected.fa")
    name_corrector.fix_name("./" + blast_folder + "/TriTrypDB-66_TbruceiLister427_2018_unknown.fa",
                            "./" + blast_folder + "/TriTrypDB-66_TbruceiLister427_2018_unknown_nameCorrected.fa")
    name_corrector.fix_name("TriTrypDB-66_TbruceiLister427_2018_VSGs.fa",
                            "TriTrypDB-66_TbruceiLister427_2018_VSGs_nameCorrected.fa")

    os.system("makeblastdb -dbtype \"nucl\" -in ./" + blast_folder +
              "/all_expressed_genomic_VSGs_nameCorrected.fa -title known_VSGs > ./" + blast_folder +
              "/blastdb_output.txt")

    os.system("blastn -query ./" + blast_folder + "/TriTrypDB-66_TbruceiLister427_2018_unknown_nameCorrected.fa" +
              " -task \"blastn\" -db ./" + blast_folder + "/all_expressed_genomic_VSGs_nameCorrected.fa" +
              " -out ./" + blast_folder + "/blastResults.tsv -outfmt \"6 qseqid sseqid pident nident length " +
              "mismatch gapopen qstart qend qlen sstart send slen evalue bitscore\" -max_target_seqs 1 -max_hsps 1" +
              " > ./" + blast_folder + "/blast_output.txt")
    print("Blast finished")


def main():
    """Execute the functions"""
    # obtain_blast_tsv()
    # put this through the blast_result_filter.R
    # returns filtered blast file where 80% of the query length is covered and the alignment has
    # a bitscore of at least 500

    # combine all Lister427 VSGs into one FASTA
    genome_prep.gff_filter("./TriTrypDB-66_TbruceiLister427_2018_unknown.gff",
                           "./" + blast_folder + "/blastResults_filtered.tsv",
                           "./TriTrypDB-66_TbruceiLister427_2018_unknown_filtered.gff")
    vsg_sequence_pull.gff_sequence_puller("TriTrypDB-66_TbruceiLister427_2018_unknown_filtered.gff",
                                          "./genome_2018_associatedFiles/" +
                                          "TriTrypDB-66_TbruceiLister427_2018_Genome.fasta",
                                          "TriTrypDB-66_TbruceiLister427_2018_unknown_filtered.fa")
    name_corrector.fix_name("TriTrypDB-66_TbruceiLister427_2018_unknown_filtered.fa",
                            "TriTrypDB-66_TbruceiLister427_2018_unknown_filtered_nameCorrected.fa")

    os.system("cat TriTrypDB-66_TbruceiLister427_2018_VSGs_nameCorrected.fa " +
              "TriTrypDB-66_TbruceiLister427_2018_unknown_filtered_nameCorrected.fa " +
              "./cross_VSGnomes/vsgs_tb427*.fa > all_Lister427VSGs.fa")
    # remove duplicates with cd-hit
    # run this on a terminal with cd-hit loaded into the path
    # os.system("cd-hit-est -i all_Lister427VSGs.fa -o all_unique_Lister427VSGs.fa -c 1.0 -n 8 -M 16000 -d 0")

    # check that all VSGs in cluster file end up picking the 2018 VSGs before duplicate are added back
    # because 2018 VSGs have a position and we know they are actually in the genome twice, we are only
    # adding back VSGs with evidence for duplication
    genome_prep.check_2018_selected_as_master("all_unique_Lister427VSGs.fa.clstr")

    # make list of all duplicates found within genome. will check them later for diversification
    # within alt family members
    # run with both files
    genome_prep.add_duplicates_back("all_unique_Lister427VSGs.fa",
                                    "TriTrypDB-66_TbruceiLister427_2018_VSGs_nameCorrected.fa",
                                    "TriTrypDB-66_TbruceiLister427_2018_VSGs_nameCorrected_duplicates.fa")
    genome_prep.add_duplicates_back("all_unique_Lister427VSGs.fa",
                                    "TriTrypDB-66_TbruceiLister427_2018_unknown_filtered_nameCorrected.fa",
                                    "TriTrypDB-66_TbruceiLister427_2018_unknown_filtered_nameCorrected_duplicates.fa")

    # make the final fasta file with duplicates added back
    # all Lister427 VSGs - many duplicates included, but are filtered out here.
    # add back VSGs which are found within 2018 genome at distinct positions but are identical
    os.system("cat all_unique_Lister427VSGs.fa TriTrypDB-66_TbruceiLister427_2018_VSGs_nameCorrected_duplicates.fa " +
              "TriTrypDB-66_TbruceiLister427_2018_unknown_filtered_nameCorrected_duplicates.fa > " +
              "all_unique_posDuplicate_Lister427VSGs.fa")
    with open("all_unique_posDuplicate_Lister427VSGs.fa", "r") as input, \
         open("all_unique_posDuplicate_Lister427VSGsonly.fa", "w") as output:
        skip = False
        for line in input:
            if ">Tb427VSG" in line and "unitig" not in line:
                skip = True
            elif skip:
                skip = False
            else:
                output.write(line)

    # # prep for genome
    # # use script genome_plotting.R to draw VSGs onto megabase chrs
    # genome_prep.genome_coords("./genome_2018_associatedFiles/TriTrypDB-66_TbruceiLister427_2018_Genome.fasta",
    #                           "genome_2018_coords.tsv")




    print("finished script")


if __name__ == '__main__':
    main()
