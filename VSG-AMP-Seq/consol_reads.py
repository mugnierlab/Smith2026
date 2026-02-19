"""VSG-AMP-Seq: This set of functions is for consolidating reads using cd-hit
   These commands rely on the os and are written for linux
   Smith 2024
   Author: Jaclyn Smith
"""
import glob
import os
import subprocess
import csv
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import consol_supplement


def umi_extract(file_name, out_put_folder):
    """
    This function isolates the umis from the titles of the reads and puts them into a
    FASTA file
    Args:
        file_name: a string, the read file to extract the UMIs from
        out_put_folder: a string, the output file name

    """
    file_name_reds = file_name.split(".gz")[0]
    out = out_put_folder + file_name.strip().split("/")[2].split("_1.")[0] + "_umis.fa"
    with open(file_name_reds, "r") as input_file, open(out, "w") as output:
        for (title, sequence, quality) in FastqGeneralIterator(input_file):
            new_title = title.replace(" ", "")
            output.write(">" + new_title + "\n")
            output.write(title.split(":")[-1].split("\\1")[0] + "\n")


def file_splitter(file_name, split_num):
    """
    Uses the os to split the UMI FASTA into multiple subfiles
    Args:
        file_name: a string, the name of the file to be split
        split_num: an int, the number of lines to put into each file
    """
    os.system("split -l " + str(split_num) + "  " + file_name + " ./split_files/")
    os.system("mv ./split_files/smaller.fa ./split_files/tmp/")
    for entry in glob.glob("./split_files/a*"):
        print(entry)
        os.rename(entry, entry + ".fa")
        os.system("cd-hit-est -i " + entry + ".fa -o " + entry + " -c 1.0 -n 8 -M 16000 -d 0")

    os.system("mv ./split_files/*.fa ./split_files/tmp/")
    os.system("mv ./split_files/*.clstr ./split_files/tmp/")
    os.system("cat split_files/* > ./split_files/smaller.fa")
    os.system("mv ./split_files/a* ./split_files/tmp/")


def file_splitter_3R(file_name):
    """
    Uses the os to split the UMI FASTA into multiple subfiles
    specifically set up for larger 3R files to regroup subfiles into different groups
    making it more likely that the same UMIs are in the same file and can be removed
    without reordering the UMI list
    Args:
        file_name: a string, the name of the file to be split
        split_num: an int, the number of lines to put into each file
    """
    os.system("split -l 4000000 " + file_name + " ./split_files/")
    os.system("mv " + file_name + " ./split_files/tmp/")
    prev_file = ""
    i = 0
    my_list = glob.glob("./split_files/a*")
    first_file = my_list.pop(0)
    for entry in my_list:
        print(entry)
        os.rename(entry, entry + ".fa")
        if i % 2 == 0:
            prev_file = entry + ".fa"
        elif i % 2 == 1:
            os.system("cat " + prev_file + " " + entry + ".fa > " + prev_file.split(".fa")[0] + "_" +\
                      entry.split("/")[2] + ".fa")
            os.system("mv " + prev_file + " ./split_files/tmp/")
            os.system("mv " + entry + ".fa ./split_files/tmp/")
        i += 1

    for reentry in glob.glob("./split_files/a*.fa"):
        os.system("cd-hit-est -i " + reentry + " -o " + reentry.split(".fa")[0] + " -c 1.0 -n 8 -M 16000 -d 0")
    os.system("mv ./split_files/*.fa ./split_files/tmp/")
    os.system("mv ./split_files/*.clstr ./split_files/tmp/")
    os.system("cat split_files/* > ./split_files/smaller.fa")
    os.system("mv ./split_files/a* ./split_files/tmp/")


def file_breakout(file_name):
    """
    Group reads by UMI, reduces UMIs to 92% similar groups, obtains consensus sequence
    from all reads. these are broken into smaller files and groups such that files are managable with
    little memory and within the limit of cd-hit
    Args:
        file_name: a string, the name of the file to be consolidated

    """
    info = file_name.strip().split("/")[2].split("_")
    print(info)
    sample_info = info[0] + "_" + info[1] + "_" + \
                  info[2] + "_" + info[3] + "_" + info[4]  # mouse, day, genotype, barcode, primer

    if os.path.exists("./consol_reads/" + sample_info + "_consol_1.fq"):
        return ()

    file_splitter(file_name, 2000000)

    count_check = subprocess.check_output("wc -l ./split_files/smaller.fa", shell=True)
    count = int(str(count_check).split(" ")[0].split("'")[1])
    if count / 8000000 >= 1:
        file_splitter("./split_files/smaller.fa", 8000000)
        count_check = subprocess.check_output("wc -l ./split_files/smaller.fa", shell=True)
        count = int(str(count_check).split(" ")[0].split("'")[1])
        if count / 8000000 >= 1:
            print("testHere... this was the last file you were on: " + file_name)
            file_splitter_3R("./split_files/smaller.fa")
            count_check = subprocess.check_output("wc -l ./split_files/smaller.fa", shell=True)
            count = int(str(count_check).split(" ")[0].split("'")[1])
        if count / 8000000 >= 1:
            file_splitter("./split_files/smaller.fa", 12000000)
            count_check = subprocess.check_output("wc -l ./split_files/smaller.fa", shell=True)
            count = int(str(count_check).split(" ")[0].split("'")[1])
        if count / 8000000 >= 1:
            file_splitter("./split_files/smaller.fa", 16000000)

    os.system("cd-hit-est -i split_files/smaller.fa -o split_files/master.fa -c 0.92 -n 8 -M 16000 -d 0")
    consensus_umis = {}
    with open("split_files/master.fa", "r") as umis:
        for line in umis:
            if ">" in line:
                consensus_umis[line.strip().split(":")[-1].split("\\1")[0]] = 0
    with open("split_files/master.fa.clstr", "r") as clusters:
        umi_codes = {}
        umi_group = []
        for line in clusters:
            if ">Cluster" in line:
                for item in umi_group:
                    umi_codes[item] = leader
                umi_group = []
                leader = ""
            else:
                if "*" in line:
                    leader = line.strip().split(":")[-1].split("\\1")[0]
                    umi_group.append(leader)
                else:
                    umi_group.append(line.strip().split(":")[-1].split("\\1")[0])
        for item in umi_group:
            umi_codes[item] = leader

    total = 0
    reads_incorporated = 0

    with open(file_name, "r") as file:
        for line in file:
            if ">" in line:
                try:
                    consensus_umis[line.split(":")[-1].split("\\1")[0]] += 1
                except KeyError:
                    try:
                        coded_umi = umi_codes[line.split(":")[-1].split("\\1")[0]]
                    except KeyError:
                        coded_umi = umi_codes[consol_supplement.rev_complement(line.split(":")[-1].split("\\1")[0])]
                    consensus_umis[coded_umi] += 1
            total += 1
    multiple_umis = []
    for item in consensus_umis:
        if consensus_umis[item] > 2:
            multiple_umis.append(item)
    my_dict = {}

    i = 1
    while i <= len(multiple_umis):
        my_dict[multiple_umis[i - 1]] = []
        if i % 20000 == 0:
            with open("./demultiplexed_reads/" + sample_info + "_1.fq", "r") as read1, open(
                    "./demultiplexed_reads/" + sample_info + "_2.fq", "r") as read2:
                for (title_1, sequence_1, quality_1), (title_2, sequence_2, quality_2) in zip(
                        FastqGeneralIterator(read1),
                        FastqGeneralIterator(read2)):
                    read_umi = title_1.strip().split(":")[-1].split("\\1")[0]
                    try:
                        final_umi = umi_codes[read_umi]
                    except KeyError:
                        final_umi = umi_codes[consol_supplement.rev_complement(read_umi)]
                    try:
                        title = title_1 + "_" + title_2
                        sequence = sequence_1 + "_" + sequence_2
                        quality = quality_1 + "_" + quality_2
                        my_dict[final_umi].append((title, sequence, quality))
                        reads_incorporated += 1
                    except KeyError:
                        pass

                with open("./consol_reads/" + sample_info + "_consol_1.fq", "a") as read1_out, \
                        open("./consol_reads/" + sample_info + "_consol_2.fq", "a") as read2_out:  # output files
                    for item in my_dict:
                        consol_read, name_record,\
                        bad_seq1_count, bad_seq2_count = consol_supplement.consensus_sequence(my_dict[item], 2, 10)
                        title_1 = consol_read[0] + "_1"
                        title_2 = consol_read[0] + "_2"
                        sequence_1 = consol_read[1].split("_")[0]
                        sequence_2 = consol_read[1].split("_")[1]
                        quality_1 = consol_read[2].split("_")[0]
                        quality_2 = consol_read[2].split("_")[1]
                        read1_out.write("@%s\n%s\n+\n%s\n" % (title_1, sequence_1, quality_1))
                        read2_out.write("@%s\n%s\n+\n%s\n" % (title_2, sequence_2, quality_2))
            my_dict = {}
        i += 1

    with open("./demultiplexed_reads/" + sample_info + "_1.fq", "r") as read1, open(
            "./demultiplexed_reads/" + sample_info + "_2.fq", "r") as read2:
        for (title_1, sequence_1, quality_1), (title_2, sequence_2, quality_2) in zip(
                FastqGeneralIterator(read1),
                FastqGeneralIterator(read2)):
            read_umi = title_1.strip().split(":")[-1].split("\\1")[0]
            try:
                final_umi = umi_codes[read_umi]
            except KeyError:
                final_umi = umi_codes[consol_supplement.rev_complement(read_umi)]
            try:
                title = title_1 + "_" + title_2
                sequence = sequence_1 + "_" + sequence_2
                quality = quality_1 + "_" + quality_2
                my_dict[final_umi].append((title, sequence, quality))
                reads_incorporated += 1
            except KeyError:
                pass

        with open("./consol_reads/" + sample_info + "_consol_1.fq", "a") as read1_out, \
                open("./consol_reads/" + sample_info + "_consol_2.fq", "a") as read2_out:  # output files
            for item in my_dict:
                consol_read, name_record,\
                bad_seq1_count, bad_seq2_count = consol_supplement.consensus_sequence(my_dict[item], 2, 10)
                title_1 = consol_read[0] + "_1"
                title_2 = consol_read[0] + "_2"
                sequence_1 = consol_read[1].split("_")[0]
                sequence_2 = consol_read[1].split("_")[1]
                quality_1 = consol_read[2].split("_")[0]
                quality_2 = consol_read[2].split("_")[1]
                read1_out.write("@%s\n%s\n+\n%s\n" % (title_1, sequence_1, quality_1))
                read2_out.write("@%s\n%s\n+\n%s\n" % (title_2, sequence_2, quality_2))

    with open("consolidation_stats.txt", "a") as statistics:
        writer = csv.writer(statistics)
        writer.writerow([info[0], info[1], info[2], info[4], total, len(multiple_umis), reads_incorporated,
                         total - reads_incorporated])


def main():
    """Execute the functions"""
    folder_umi_fastas = "split_files"
    try:
        os.mkdir(folder_umi_fastas)
        print("folder '{}' created ".format(folder_umi_fastas))
    except FileExistsError:
        print("folder {} already exists".format(folder_umi_fastas))

    folder_umi_fastas = "./split_files/tmp"
    try:
        os.mkdir(folder_umi_fastas)
        print("folder '{}' created ".format(folder_umi_fastas))
    except FileExistsError:
        print("folder {} already exists".format(folder_umi_fastas))

    folder_consol = "consol_reads"
    try:
        os.mkdir(folder_consol)
        print("folder '{}' created ".format(folder_consol))
    except FileExistsError:
        print("folder {} already exists".format(folder_consol))

    header = ["mouse", "day", "genotype", "primer", "total_reads",
              "num_consol_reads", "reads_consolidated", "reads_eliminated"]
    with open("consolidation_stats.txt", "a") as statistics:
        writer = csv.writer(statistics)
        writer.writerow(header)

    folder_umi_fastas = "./umi_fastas/"
    try:
        os.mkdir(folder_umi_fastas)
        print("folder '{}' created ".format(folder_umi_fastas))
    except FileExistsError:
        print("folder {} already exists".format(folder_umi_fastas))

    for file in glob.glob("demultiplexed_reads/*_1.fq.gz"):

        folder = "split_files"
        try:
            os.mkdir(folder)
            print("folder '{}' created ".format(folder))
        except FileExistsError:
            print("folder {} already exists".format(folder))

        folder = "./split_files/tmp"
        try:
            os.mkdir(folder)
            print("folder '{}' created ".format(folder))
        except FileExistsError:
            print("folder {} already exists".format(folder))

        file_info = file.split("_1.fq.gz")[0]
        os.system("gunzip ./" + file_info + "*.gz")
        file_umis = file.split("/")[1].split("_1.")[0] + "_umis.fa"
        umi_extract("./" + file_info + "_1.fq", folder_umi_fastas)
        print("./umi_fastas/" + file_umis)
        file_breakout("./umi_fastas/" + file_umis)
        os.system("rm -r split_files")
        os.system("gzip ./" + file_info + "*.fq")


if __name__ == '__main__':
    main()
