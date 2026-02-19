"""VSG-AMP-Seq: additional functions
   Smith 2024
   Author: Jaclyn Smith
"""

gen_code = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N", "X": "X", "c": "g", "g": "c",
            "a": "t", "t": "a", "n": "n"}

protein_code = {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
                "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
                "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
                "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
                "TAT": "Y", "TAC": "Y", "TAA": "1", "TAG": "1", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
                "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
                "TGT": "C", "TGC": "C", "TGA": "1", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
                "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
                "TCN": "S", "CTN": "L", "GTN": "V", "CCN": "P", "ACN": "T", "GCN": "A", "CGN": "R", "GGN": "G"}


def read_in_primers(primers):
    """Primer Read in

    This reads in primers from file into a dictionary; also creates dictionary for
    vsg predicted sequences

    Args:
        primers: tab delimited file with primer names followed by sequences followed by
        the predicted 150 bp sequence from the VSG of interest
    """
    # dict with primers
    primer_dictionary = {}
    # dict with predicted seq
    pred_seq_dict = {}
    with open(primers, "r") as primer_file:
        for entry in primer_file:
            primer = entry.strip().split("\t")
            # key = sequence, value = primer name
            primer_dictionary[primer[1]] = primer[0]
            pred_seq_dict[primer[0]] = primer[2]
    return primer_dictionary, pred_seq_dict


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
        new_base = gen_code[position]
        comp_read = comp_read + new_base
    rev_comp_read = comp_read[::-1]
    return rev_comp_read


def main():
    """Execute the functions"""
    print("none")


if __name__ == '__main__':
    main()
