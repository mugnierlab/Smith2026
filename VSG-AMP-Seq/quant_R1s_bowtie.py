"""VSG-AMP-Seq: quantification of the R1 unanchored reads
   Smith 2024
   Author: Jaclyn Smith
"""
import os
import glob

locations = {"141": 243,
             "300R": 369,
             "592": 694,
             "792": 894,
             "909R": 978,
             "1357": 1459}


def run_bowtie(file_name, path):
    antat_folder = "./bowtie_files/"
    try:
        os.mkdir(antat_folder)
        print("folder '{}' created ".format(antat_folder))
    except FileExistsError:
        print("folder {} already exists".format(antat_folder))

    output_folder = "./sam_alignments/"
    try:
        os.mkdir(output_folder)
        print("folder '{}' created ".format(output_folder))
    except FileExistsError:
        print("folder {} already exists".format(output_folder))

    if not os.path.exists(antat_folder + "antat" + ".1.ebwt"):
        os.system("bowtie-build " + file_name + " " + antat_folder + "antat")

    for file in glob.glob(path + "*.fq.gz"):
        os.system("gunzip " + file)
        sample_info = file.split(".fq")[0].split("/")[-1]
        os.system("bowtie -v 2 -S --no-unal " + antat_folder + "antat -q " + file.split(".gz")[0] + " &> " +
                  output_folder + sample_info + ".sam")
        os.system("gzip " + file.split(".gz")[0])


def sam_location():
    for item in glob.glob("./sam_alignments/*.sam"):
        guide = item.split("_")[3].split(".")[1]
        position = locations[guide]
        space = 250
        start_range = position - space
        end_range = position + space
        counted_reads = 0
        with open(item, "r") as sam_file, open("output", "a") as stats:
            for line in sam_file:
                if line.startswith("@"):
                    pass
                elif line.startswith("#"):
                    pass
                elif line.startswith("Reported"):
                    pass
                elif line.startswith("Setting"):
                    pass
                else:
                    len_read = len(line.split("\t")[9])
                    if start_range <= int(line.split("\t")[3]) <= end_range:
                        counted_reads += 1
                    elif start_range <= int(line.split("\t")[3]) + len_read <= end_range:
                        counted_reads += 1
            stats.write(item.split(".sam")[0] + "\t" + str(counted_reads) + "\n")


def main():
    antat_fasta = "antat.fasta"
    run_bowtie(antat_fasta, "../single_marker_228Mechanism/")
    sam_location()


if __name__ == '__main__':
    main()
