"""VSG-AMP-Seq: This is a group of functions to supplement the vsg_align module
   Smith 2024
   Author: Jaclyn Smith"""
import glob
import global_target as gbl


def read_in_candidateVSGs(primer, status):
    """
    Reads in sorted donor VSGs according to primer presence
    Args:
        primer: a string, the primer identity
        status: a string, either "contains" or "absent"
    Returns:
        vsg_candidate_dict: a dictionary
            key = name of the VSG
            value = sequence
    """
    # look for reads aligning to alternative VSGs with primer sequence
    # Can put in any FASTA into the folder donor_fastas for VSG readin
    folder_donors = "./donor_fastas/"
    vsg_candidate_dict = {}
    for vsg_file in glob.glob(folder_donors + "*" + primer + "_" + status + ".fa"):
        with open(vsg_file, "r") as subset_VSGs:
            title = ""
            sequence = ""
            for line in subset_VSGs:
                if ">" in line:
                    title = line.strip().split(">")[1].split(" :")[0]
                else:
                    sequence = line
                if title != "" and sequence != "":
                    vsg_candidate_dict[title] = sequence

                    title = ""
                    sequence = ""
    return vsg_candidate_dict


def rev_complement(read):
    """
    Creates the reverse complement of a sequence
    Args:
        read: a string, a DNA sequence

    Returns:
        rev_comp_read: a string, the reverse complement of read
    """
    comp_read = ""
    for position in read:
        new_base = gbl.gen_code[position]
        comp_read = comp_read + new_base
    rev_comp_read = comp_read[::-1]
    return rev_comp_read


def lowercase_primer(primer, sequence_2):
    """
    Makes the primer portion of the sequence lowercase
    Args:
        primer: string, the primer identity
        sequence_2: string, the read sequence for R2
    Returns:
        sequence_2: string, the updated string with the primer lowercase
    """
    primer_seq = gbl.rev_primer_dict[primer]
    primer_lower = len(primer_seq)
    sequence_2 = sequence_2[:primer_lower].lower() + sequence_2[primer_lower:]
    return sequence_2
