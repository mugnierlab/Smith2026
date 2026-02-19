"""VSG-AMP-Seq: This set of functions demultiplexes reads by sample with low memory usage
   Smith 2024
   Author: Jaclyn Smith
"""
import os
import glob
from contextlib import ExitStack
import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import Levenshtein as lev


def barcodes(barcode_file):
    """Barcode Read in

    This reads in barcodes from file into a dictionary

    Args:
        barcode_file: tab-delimited .txt file with sample name,
        i7 barcode sequence, and i5 barcode sequence
    """
    barcode_dict = {}

    # barcode is a tab delimited file with names followed by sequence
    with open(barcode_file, "r") as bc_file:
        for entry in bc_file:
            barcode = entry.strip().split("\t")
            # key = i7+i5, value = sample name
            barcode_dict[barcode[1] + "+" + barcode[2]] = barcode[0]
    return barcode_dict


def barcode_errors(barcode_dict):
    """
    Barcode error Calculation
    Compare used barcode to determine the allowable mismatch to still
    identify a unique barcode.
    Args:
        barcode_dict: dictionary of barcodes

    Returns:
        reasonable_error: int, error that should be allowed

    """
    # initialize variables
    lowest_dist = 16
    extra_space = 2

    # loop through barcode_dict keys twice
    for entry in barcode_dict:
        for entry2 in barcode_dict:
            # identify differences between barcodes
            dist = (lev.distance(entry, entry2))
            # identify most similar barcodes
            if dist != 0:
                lowest_dist = min(dist, lowest_dist)
    # calculate a reasonable error rate to allow
    reasonable_error = lowest_dist - extra_space
    return reasonable_error


def demultiplex(barcode_dict, error):
    """
    Demultiplex by experiments

    This sorts the reads into files based on bar codes present.
    Each experiment will have one file per primer.

    Args:
        barcode_dict: dictionary from barcodes() with i7+i5 as key
        and sample names as values
        error: int, allowed mismatches true barcode
    """
    # initialize variables

    total_unknown_barcodes = 0
    total_still_unsorted = 0
    unused_barcodes = {}
    UMI_len = 25
    data = []
    data_by_primer = []

    folder_output = "demultiplexed_reads"
    try:
        os.mkdir(folder_output)
        print("folder '{}' created ".format(folder_output))
    except FileExistsError:
        print("folder {} already exists".format(folder_output))

    # loop through all files in folder trimmed_reads
    for trimm_file in glob.glob("./removed_spacer_reads/*_1.fq"):
    #for trimm_file in glob.glob("./consol_reads/*_1.fq"):
        r2_trimm_file = trimm_file.split("_1.fq")[0] + "_2.fq"

        # read in sequences from file into dictionary as tuples
        with open(trimm_file, "r") as seq_file1, open(r2_trimm_file, "r") as seq_file2:
            # keys = barcodes from barcode_dict, values = empty list, for reads tuples
            multiplex_counts = {key: 0 for key in barcode_dict}
            primer = trimm_file.split("/")[2].split("_")[1]
            barcode_positions = {}
            output_files = []
            counter = 0
            folder = "./demultiplexed_reads/"
            for barcode in barcode_dict:
                barcode_positions[barcode] = counter
                counter += 2
                file_name_r1 = folder + barcode_dict[barcode] + "_" + barcode + "_" + primer + "_1.fq"
                file_name_r2 = folder + barcode_dict[barcode] + "_" + barcode + "_" + primer + "_2.fq"

                output_files.append(file_name_r1)
                output_files.append(file_name_r2)
            output_files.append(folder + "unassigned_" + primer + "_1.fq")
            output_files.append(folder + "unassigned_" + primer + "_2.fq")

            with ExitStack() as stack:
                files = [stack.enter_context(open(fname, "w")) for fname in output_files]
                for (title_1, sequence_1, quality_1), (title_2, sequence_2, quality_2) in zip(
                        FastqGeneralIterator(seq_file1),
                        FastqGeneralIterator(seq_file2)):
                    barcode = title_1.split(":")[-1][:-UMI_len]
                    unique_molecular_index = title_1[-UMI_len:]
                    try:
                        multiplex_counts[barcode] += 1

                        # alter read name, add unique molecular index instead of
                        # barcode, add \ read num on end of read name.
                        mod_read_name_r1 = title_1.rsplit(":", 1)[0] + ":" + \
                                        unique_molecular_index + "\\" + \
                                        title_1.split(" ")[1].split(":")[0]
                        mod_read_name_r2 = title_2.rsplit(":", 1)[0] + ":" + \
                                        unique_molecular_index + "\\" + \
                                        title_2.split(" ")[1].split(":")[0]
                        # proper fastq format
                        files[barcode_positions[barcode]].write(
                            "@%s\n%s\n+\n%s\n" % (mod_read_name_r1, sequence_1, quality_1))
                        files[barcode_positions[barcode]+1].write(
                            "@%s\n%s\n+\n%s\n" % (mod_read_name_r2, sequence_2, quality_2))

                    # keyError here is for reads from previous runs/run errors
                    # some barcodes here are not in the barcode file, which will throw an error
                    except KeyError:
                        # count the reads with unknown barcodes
                        total_unknown_barcodes = total_unknown_barcodes + 1
                        no_id = True
                        for entry in multiplex_counts:
                            dist = lev.distance(entry, barcode)
                            if dist <= error:
                                multiplex_counts[entry] += 1

                                # alter read name, add unique molecular index instead of
                                # barcode, add \ read num on end of read name.
                                mod_read_name_r1 = title_1.rsplit(":", 1)[0] + ":" + \
                                                   unique_molecular_index + "\\" + \
                                                   title_1.split(" ")[1].split(":")[0]
                                mod_read_name_r2 = title_2.rsplit(":", 1)[0] + ":" + \
                                                   unique_molecular_index + "\\" + \
                                                   title_2.split(" ")[1].split(":")[0]
                                # proper fastq format
                                files[barcode_positions[entry]].write(
                                    "@%s\n%s\n+\n%s\n" % (mod_read_name_r1, sequence_1, quality_1))
                                files[barcode_positions[entry]+1].write(
                                    "@%s\n%s\n+\n%s\n" % (mod_read_name_r2, sequence_2, quality_2))
                                no_id = False
                                break
                        if no_id:
                            # keep track of what these unused bar codes are and how many reads with each
                            try:
                                total_still_unsorted += 1
                                # previously seen unused barcode is counted
                                unused_barcodes[barcode] = unused_barcodes[barcode] + 1
                            except KeyError:
                                # if barcode not seen before, added to dictionary
                                unused_barcodes[barcode] = 1
                            files[len(files) - 2].write("@%s\n%s\n+\n%s\n" % (title_1, sequence_1, quality_1))
                            files[len(files) - 1].write("@%s\n%s\n+\n%s\n" % (title_2, sequence_2, quality_2))

            for demulti_file in multiplex_counts:
                info = barcode_dict[demulti_file].split("_")
                mouse = info[0]
                day = info[1]
                genotype = ""
                direction = primer[1]
                data.append((mouse, day, genotype, direction, multiplex_counts[demulti_file]))
                data_by_primer.append((mouse, day, genotype, direction, multiplex_counts[demulti_file]))
    print(str(total_unknown_barcodes) + " reads had unknown bar codes.")
    print(str(total_still_unsorted) + " reads still cannot be sorted")
    print("The following unknown barcode was observed the most:")
    most_unused_barcode = max(unused_barcodes, key=unused_barcodes.get)
    print(most_unused_barcode + " " + str(unused_barcodes[most_unused_barcode]))
    print("Total unknown barcodes: " + str(len(unused_barcodes)))

    summary_data = pd.DataFrame(data, columns=["mouse", "day", "genotype", "direction", "count"])
    total_reads = summary_data.sum()["count"]
    ratio_per_sample = summary_data.groupby(["mouse", "day", "direction"]).sum()
    sample_ratios = ratio_per_sample.pivot_table(index=["mouse", "day"], columns="direction", values="count")
    sample_ratios["total_reads_sample"] = sample_ratios["F"] + sample_ratios["R"]
    sample_ratios["percentF"] = sample_ratios["F"] / sample_ratios["total_reads_sample"] * 100
    sample_ratios["percentR"] = sample_ratios["R"] / sample_ratios["total_reads_sample"] * 100
    print("Summary Data:")
    print("This is summary data for each sample F vs R")
    print(sample_ratios)
    sample_ratios.to_csv("sample_ratios.csv")

    # summarise data by sample, results printed and written to CSV
    total_reads_per_sample = ratio_per_sample.groupby(["mouse", "day"]).sum()
    total_reads_per_sample["sample_ratio"] = total_reads_per_sample["count"] / total_reads * 100
    print("")
    print("Num reads per sample:")
    print("This is summary data for each sample, post primer sort.")
    print(total_reads_per_sample)
    total_reads_per_sample.to_csv("total_reads_per_sample.csv")

    # summarise data by primer, results written to CSV
    primer_summary = pd.DataFrame(data_by_primer, columns=["mouse", "day", "genotype", "primer", "count"])
    primer_summary.to_csv("readsByPrimer.csv")


def main():
    """Execute the functions"""

    bc = barcodes("JAC22-6b_barcodes.txt")
    allowed_error = barcode_errors(bc)
    demultiplex(bc, allowed_error)


if __name__ == '__main__':
    main()
