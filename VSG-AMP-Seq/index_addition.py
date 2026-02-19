"""VSG-AMP-seq: Adds index sequence to titles of sequences
   Smith 2024
   Author: Jaclyn Smith"""

import os
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import global_target as gbl


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


def index_add(read1, read2, index1, index2, sequencer="MiSeq"):
    """
    Add the indexes onto the read names

    This adds indexes from the index files onto the read files

    Args:
        read1: total fastq file for read1
        read2: total fastq file for read2
        index1: index fastq file contains i7
        index2: index fastq file contains i5 + UMI and sample names as values
        sequencer: a string, the type of machine the data was generated with.
            I2 may need to be inverted depending on sequencing chemistry.
    """

    folder_reads = "index_reads"
    try:
        os.mkdir(folder_reads)
        print("folder '{}' created ".format(folder_reads))
    except FileExistsError:
        print("folder {} already exists".format(folder_reads))

    folder = "./index_reads/"
    new_output_1 = folder + "index_R1.fastq"
    new_output_2 = folder + "index_R2.fastq"

    # open all files. Loop through reads and indexes simultaneously.
    with open(read1) as r1, open(read2) as r2, open(index1) as i1, open(index2) as i2, \
            open(new_output_1, "w") as out1, open(new_output_2, "w") as out2:
        for (title_1, sequence_1, quality_1), \
            (title_2, sequence_2, quality_2), \
            (title_i1, sequence_i1, quality_i1), \
            (title_i2, sequence_i2, quality_i2) in zip(
            FastqGeneralIterator(r1),
            FastqGeneralIterator(r2),
            FastqGeneralIterator(i1),
            FastqGeneralIterator(i2)):

            if sequencer == "NovaSeq":
                update_i2 = rev_complement(sequence_i2)
            elif sequencer == "MiSeq":
                update_i2 = sequence_i2

            # alter the read names to include index
            mod_read_name_1 = title_1.rsplit(":", 1)[0] + ":" + \
                              sequence_i1 + "+" + update_i2
            mod_read_name_2 = title_2.rsplit(":", 1)[0] + ":" + \
                              sequence_i1 + "+" + update_i2

            # write to output
            out1.write("@%s\n%s\n+\n%s\n" % (mod_read_name_1, sequence_1, quality_1))
            out2.write("@%s\n%s\n+\n%s\n" % (mod_read_name_2, sequence_2, quality_2))


def main():
    """Execute the functions"""
    rev_complement("AGGTTA")
    index_add("tissues_R1.fastq", "tissues_R2.fastq",
              "tissues_I1.fastq", "tissues_I2.fastq", "NovaSeq")


if __name__ == '__main__':
    main()
