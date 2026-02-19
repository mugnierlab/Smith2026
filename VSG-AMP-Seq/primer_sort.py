"""VSG-AMP-Seq: This function sorts reads by primer and reports sorting results.
   Smith 2024
   Author: Jaclyn Smith"""

import os

gen_code = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N", "X": "X", "c": "G", "g": "C",
            "a": "T", "t": "A"}


def rev_complement(read):
    """
    Reverse complement function - for DNA seqs
    Args:
        read: string, DNA sequence

    Returns:
        rev_comp_read: string, reverse and complemented DNA seq

    """
    comp_read = ""
    for position in read:
        new_base = gen_code[position]
        comp_read = comp_read + new_base
    rev_comp_read = comp_read[::-1]
    return rev_comp_read


def sort_sequences(expt_name, primer_fasta, seq_file_1, seq_file_2):
    """Sorts sequences

    Locates primers on 5' end of read 2, sorts all paired sequences into
    destination files based on primer presence

    Args:
        expt_name: experiment header for output file name
        primer_fasta: fasta file containing primers for target
        seq_file_1: First read fastq file
        seq_file_2: Second read fastq file
    """
    # initialize flags/variables

    # universal cutadapt command variables
    # overlap 10, for spacer length of 7 + default length 3
    # yields properly sorted reads
    # here seq2 is purposely listed first!!!! read used for demultiplexing (R2) must be first, this is consistent for
    # primer flags, input and output
    cut_flags = "--action=retain --overlap 10 "
    input_files = seq_file_2 + " " + seq_file_1

    # create folder for all sorted read files
    folder_sorted_reads = "sorted_reads"
    try:
        os.mkdir(folder_sorted_reads)
        print("folder '{}' created ".format(folder_sorted_reads))
    except FileExistsError:
        print("folder {} already exists".format(folder_sorted_reads))
    folder = "./sorted_reads/"

    # creating variables to put into cutadapt command
    # search for sequence on 5' end of read2, but R2 is first, so lowercase letter is used!!
    # g for 5' end of read, no ^ because read is unanchored (after spacer was added)
    # primers are supplied in a fasta format; names of primers are fasta headers
    # this is important because all primer sequences are compared simultaneously and best match is chosen
    # for input and output, R2 is first bc that contains the primer used to demutliplex
    # and to deal with this, uppercase letter changed to lowercase bc read2 is in read1 position now
    adapter_command = "-g file:" + primer_fasta + " "
    output_read_1 = "-o " + folder + "{name}_read2.fq "
    output_read_2 = "-p " + folder + "{name}_read1.fq "
    unsorted_file_output = "--untrimmed-output orphan_read2.fq" + \
        " --untrimmed-paired-output orphan_read1.fq "

    cut_command = "cutadapt " + adapter_command + cut_flags + output_read_1 + \
                  output_read_2 + unsorted_file_output + input_files
    os.system(cut_command + "> sort_reads_report_" + expt_name + ".txt")


def main():
    """Execute the functions"""
    sort_sequences("JAC22-11a", "antat_primers.fasta", "./index_reads/index_R1.fastq", "./index_reads/index_R2.fastq")


if __name__ == '__main__':
    main()
