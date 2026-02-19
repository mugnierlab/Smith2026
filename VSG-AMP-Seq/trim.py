"""VSG-AMP-Seq: This set of functions trims reads with trimgalore and removes R1 spacer
   Smith 2024
   Author: Jaclyn Smith
"""
import os
import glob
import shutil
import aux_functions


def trim(folder_name, primer_dict, VSG_name):
    """
    Trimming reads

    This executes trim_galore on a folder containing sorted reads.
    It trim reads based on quality and remove adapters which might be
    on the 3' ends of reads to leave just the hypothetical antat sequence.
    """
    # Runs a quality and adapter based trimming process on the reads

    # checks that input folder_name is a folder
    if os.path.isdir(folder_name):
        print(folder_name + " exists.")
    else:
        print("WARNING:The folder " + folder_name + " does not exist.")

    # trim_galore command variables
    # quality trim is 1st, default settings
    trim_flags = "--dont_gzip --paired --trim1 "
    # adapters, trimmed from 3' end:
    # read1 = universal antat primer addition: ATCACCGA
    # ^ this cuts off all the invariant sequence added by each Antat primer
    # read2 = illumia primer: AGATCGGAAGAGC
    # these primers may be read through during sequencing
    read1_adapter = "ATCACCGACTGCC"
    read2_adapter = "AGATCGGAAGAGC"
    adapters = "--adapter " + read1_adapter + " --adapter2 " + read2_adapter + " "

    # create folder for trim_galore output files
    folder_output = "trim_galore"
    try:
        os.mkdir(folder_output)
        print("folder '{}' created ".format(folder_output))
    except FileExistsError:
        print("folder {} already exists".format(folder_output))

    folder_reads = "trimmed_reads"
    try:
        os.mkdir(folder_reads)
        print("folder '{}' created ".format(folder_reads))
    except FileExistsError:
        print("folder {} already exists".format(folder_reads))

    for primer_name in primer_dict.values():
        # output file name
        output_name = "--basename " + VSG_name + "_" + primer_name + " "

        # input files
        input_reads = folder_name + VSG_name + "_" + primer_name + "_read1.fq " + \
                      folder_name + VSG_name + "_" + primer_name + "_read2.fq"

        # trim_galore command
        trim_command = "trim_galore " + trim_flags + adapters + output_name + input_reads

        # run trim galore, send output to file in folder trim_galore
        print("TRIMMING")
       # must remove &> if you are on linux and not mac!!
        os.system(trim_command + " &> trim_galore/trim_output_" + primer_name + ".txt")

    # move normal trim_galore report files to trim_galore folder
    for txt_file in glob.iglob(os.path.join("./", "*_trimming_report.txt")):
        shutil.move(txt_file, "./trim_galore/")

    # move trimmed files to trimmed_reads folder
    for txt_file in glob.iglob(os.path.join("./", "*_val_*")):
        shutil.move(txt_file, "./trimmed_reads/")


def spacer_trim(primer_dictionary, VSG_name):
    """Spacer Trim

    Trims spacers from 3' end read1

    Args:
        primer_dictionary: from sort reads by primer
    """
    print("Trimming spacers from read1.")

    # initialize cutadapt folder
    folder = "cutadapt"
    try:
        os.mkdir(folder)
        print("folder '{}' created ".format(folder))
    except FileExistsError:
        print("folder {} already exists".format(folder))

    # universal cutadapt command variables
    # overlap 10, for spacer length of 7 + default length 3
    # yields properly sorted reads

    cut_flags = "--action=retain --overlap 10 "
    # I do not know if I need the overlap part, but I will leave it in case there might be an issue

    # create folder for all sorted read files
    folder_sorted_reads = "removed_spacer_reads"
    try:
        os.mkdir(folder_sorted_reads)
        print("folder '{}' created ".format(folder_sorted_reads))
    except FileExistsError:
        print("folder {} already exists".format(folder_sorted_reads))
    folder = "./trimmed_reads/"
    folder_sorted = "./removed_spacer_reads/"

    primers = {val: key for key, val in primer_dictionary.items()}

    # loop through all primers in dictionary and use as adapters for cutadapt
    for trimmed_file in glob.glob(folder + "*_1.fq"):
        print("File: " + trimmed_file)
        primer = trimmed_file.strip().split("/")[2].split("_")[1]
        seq_file_1 = folder + VSG_name + "_" + primer + "_R1_val_1.fq"
        seq_file_2 = folder + VSG_name + "_" + primer + "_R2_val_2.fq"
        input_files = seq_file_1 + " " + seq_file_2
        sequence = aux_functions.rev_complement(primers[primer])

        # creating variables to put into cutadapt command
        # search for sequence on 5' end of read2
        # lowercase letter for read1, a for 3' end of read
        adapter_command = "-a " + sequence + " "
        output_read_1 = "-o " + folder_sorted + VSG_name + "_" + primer + "_val_1.fq "
        output_read_2 = "-p " + folder_sorted + VSG_name + "_" + primer + "_val_2.fq "
        cut_command = "cutadapt " + adapter_command + cut_flags + output_read_1 + \
                      output_read_2 + input_files

        os.system(cut_command + " > cutadapt/cutadapt_output_seq1_" + VSG_name + "_" + primer + ".txt")


def main():
    """Execute the functions"""

    primer_dict, prime_seq = aux_functions.read_in_primers("antat_primers.txt")
    trim("./sorted_reads/", primer_dict, "AnTat")
    spacer_trim(primer_dict, "AnTat")


if __name__ == '__main__':
    main()
