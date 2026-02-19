"""VSG-AMP-Seq: This set of functions aligns reads to target VSG and identifies mosaic VSG reads
   Smith 2024
   Author: Jaclyn Smith"""
import os
import regex
import aux_functions
import vsg_align_supplement


def genome_primer_sort(file_name, primer_file):
    """
    Sort out potential donors that can bind the primers. Writes output to two files
        1) _contains.fa which has VSGs which bind the primer w/ 1 mismatch tolerance
        2) _absent.fa which doesn't have VSGs which bind with primer
    Args:
        file_name: name of fasta file to be sorted

    """
    # initialize variables
    title = ""
    sequence = ""
    allowed_mismatch = 1

    folder_donors = "donor_fastas"
    try:
        os.mkdir(folder_donors)
        print("folder '{}' created ".format(folder_donors))
    except FileExistsError:
        print("folder {} already exists".format(folder_donors))
    folder_donors = "./" + folder_donors + "/"

    # read in primers + sequences
    primer_dict, primer_seq_dict = aux_functions.read_in_primers(primer_file)

    for primer in primer_dict:

        if "R" in primer_dict[primer]:
            primer_seq = vsg_align_supplement.rev_complement(primer)
        else:
            primer_seq = primer

        with open(file_name, "r") as VSGenome:
            output_contains = (folder_donors + file_name.strip().split(".")[0] + "_" + primer_dict[primer] +
                               "_contains.fa")
            output_absent = (folder_donors + file_name.strip().split(".")[0] + "_" + primer_dict[primer] +
                             "_absent.fa")
            total_seqs = 0
            contains_count = 0
            absent_count = 0
            with open(output_contains, "w") as contains:
                with open(output_absent, "w") as absent:
                    for line in VSGenome:
                        if ">" in line:
                            title = line
                        else:
                            sequence = line
                            total_seqs += 1
                        if title != "" and sequence != "":
                            read_regex = ("(?b)(" + primer_seq + "){s<=" +
                                          str(allowed_mismatch) + "}")
                            primer_present = regex.search(read_regex, sequence)
                            if primer_present:
                                contains.write(title)
                                contains.write(sequence)
                                contains_count += 1

                            else:
                                absent.write(title)
                                absent.write(sequence)
                                absent_count += 1

                            title = ""
                            sequence = ""
                    print("Total number of VSGs: " + str(total_seqs))
                    print("Contains primer " + primer_dict[primer] + ": " + str(contains_count))
                    print("Primer Absent: " + str(absent_count))


def main():
    genome_primer_sort("EATRO1125_vsgs_long.fa", "antat_primers.txt")


if __name__ == '__main__':
    main()
