import os
import glob


vsg_8_seq = "GGAAAACCGAAACGATTCCAAACACTATTCGACGTTAAGCGGCAAGTGATCTATATCCCTAGGCGAGGTTAAAGGAGACACAATCGCCAAGATGCGAACAGCGAGCACAACCTTGATATTACTATTTTTGGCACTAGGATTTTTTGCCGGGACTTCACACTGCCACAAAGCGGCGCTACCAATTGCAGGGCTCACAAAAGTCTGTACTTTCAACGGCGAGTTAAAAAAGGCAGCAGCGCGCGCAACAAACCTGATCAGTAGCTACGTATCGAAACTAATCCAACTGCAGACGCTGACGGATGACATCAAACTCCTAGTAACAGAGGGACTACTAAACGCAACAGACGAAATCGCACTCATCCAAGCTGCTGGCCGTAAAGGGACAGAAGACATGATAAATACCTTGCATAACGACATGCCAAAAGCCATATACGCAGCAGCTACGTGCAGTCTTTTTGCCGGAAGGGTGGACGATTTTGTCGGCGTTTTCGCCGGAGCCAAAGGATCGGCGAACTATTGCGTCAAAGACGGCGACGGAGCTTACAACCACGCACAGGAAAACAAGCTCCAGGGGTGCCTTGAAGGCGACGGTAGTTTCCAGGGTATCTCAGAAAAGGAGCAAAACGAGCCGCCGCGACTAGGCACCAGTTATAAAGCGCTGCAAGCGGACACGGTCGGGCGCACGGCCACGGGGTCAACATGCGCCCTAACACAGCACGGCACAAACGGCGGCGAAAGCTACGCAGAGAACGGGCAAAACCCAAACGTGGTCTGGGGAAACGGGTTATTCAAAGTAGCCAACGGCGGCGCAGCAGCGGCTGTAACCGACTGGCCAATACACGCGCAGCAAACAATAACGGACTCGAAGTTTGCGGCTTGCCAAACAAAGATGCAGGAGGCGGAGAACACCCTCTTTCAGTCGTCAAACCTACAAACCAAGCTGTTGGAACTCACGAACAAGCTACCTAAAATATCAGGCAAGCTGCAGATACCAAAATCATTTTTCGACGCTAACAACCCAGCAGCCGGCTACGACGTAGATGAGGCGCAGCTTCAGCACCTACAAGCGGCGCTGGCATCAATACGAGAAACCACTGAATACAAAACCAAGAGAAAGACAAAACTTCTCACTTCTTCAAAATCACTGATATCCGCCGCAAGTAAAAAGTGTTCAGAACCACCAGCAGCGGAGACAGTTAAAAAAAGTAACACAAACAAGAAAGCAGATAAAGAATGCGATTCCAAGGGTAAAAGTGAGTGCAATGACCCCTGTGTATGGAAAGGAACCAACAACGAAGGGAAGTGCGTAGAAAAAGAGGAAGCGAAACCGGAAGAAAAGAAAGAAGAGAAATGCAAAGATAAGACAAAGGAGGGAGATTGTACTGGAAATTGCAAATGGGAGGGTGAAACTTGCAAAGATTCCTCTATTCTAGTAACCAAGAAATTCGCCCTCAGTGTGGTTTCTGCTGCATTTGTGGCCTTGCTTTTTTAAAAGAAAAATTTCCCCCCTCAATTTTTCTTGCTGAAATTTACTAAAATTCTTGCTACTTGAAAACTTCTGATATATTTTAA"
vsg_4806 = "GGAAAACCAAAACGATTCCAAACACTATTCGACGTTAAGAGGCAAGTGATCTATATCCCTAGGCGAGGTTAAAGGAGACACAATCGCCAAGATGCGAACAGCGAGCACAACCTTGATATTACTATTTTTGGCACTAGGCTTTTTTGCCGGGACTTCACACTGCCACAAAGCGGCGCTACAAATTGCAGGGCTCACAAAAGTCTGTACTTTCAACGGCGAGTTAAAAAAGGCAGCAGCGCGCAGCGCAAACCTGATCGGTATCTACGTATCGAAGCTGAGCCAACTGCAAACGCTAGAGCATAACATTAGCCTCCTAGTAATGGAAGGGGTAGTAAACGGCACCACTGATATCGCGGTCCTCAGTCACGTGGCGCGCAGAGCTACAACGGACATGGTCGATGCTTTGCTTGCAGAGATGCCAAAAGCGATATACGCAGCAGCAACGTGCAGTCTTTTTGCTGGAAGGGTGGACGATTTTGTCGGCGTTTTCGCGGCTGGGAAAGGCTCTACAAACTATTGCGTCAAAGAAGGCGGCGGCGGCTACAGCCACTCAGCCAAAAACAAACTCCAGGGGTGCCTTGAAGACGACGGTAGTTTCCAGGGTATCTCAGAAAAGGAGCAAAACGAGCCGACGCGTCTAGACAGCAGTGATAAAGCGCTGCAAGCGGACACGGCCGGCCGAACACCGACGGGGACAACATGCGTTCTAACACAGCACGGCAACAACGGTGGTGCAGCGTACGCCGAGGACGGCCAAAACCCTGACTTCCACTGGGGTAACGGCTTATTCAAGGTAGCCCACGGCGACGCAGCAGCCAAAGAAACCAACTGGCCGATGCACACGCAGGATAAACTAGCCAAAACGCAGTTTTCAGCGTGCCAAGCCAAAATGGCCGCCGTGAAAAACATCCTCCTTCAGTCAACGGCACTACAAACCAAGTTGTTGAAACTCACGACCAAACAGCCAAAAATAACAGGCAAGCTCCAAATACCAAAATCATTTTTCGACACTAACGAGCCAGCAACCAACTACGACGTGGAAGAGGCACAACTTCAGCACCTACAAGCAGCGCTGGCATCAATGCGGGAAACCACTGAATGCAAAAATAAGAGAAAAGGCAAGCTTCTCACTTCGTCAAAAACATTGATCTCCAGCGCGAGCAAAAAGTGTCAGGATTCAGCAGCAGCGGAGACAGTTAAAAAAAGTAACACAAACAAGAACGCAGATAAAGAATGCGATTCCAAGGGTAAAAGTGAGTGCAATGACCCCTGTGTATGGAAAGGAACCGACAAGGAAGGGAAGTGCGTAGAAAAAGAGGAAGCGAAACCTGAAGAAAAGAAAGAAGAGAAATGTAAAGAAAAAAAGGATGACTGCAAATCTCCGGATTGTAAGTGGGAGGGTGAAACTTGCAAAGATTACGGTTTTCTCGACAATATTAAATTTTCCTGGTGACGTCCAACACCTCATCATGGTCCGAGGGCCTGTTACTCCGTCCAGGAGGAATCCAAGTGTCCGCATTATACCATGTTCCAGAGGACTTTGGCGTGCATAAAATTTTCGGCTACCAGAGA"


def convert_bam_to_fastq(data_path):
    path_to_data = data_path + "/bam_files/bam_pass/"
    for item in glob.glob(path_to_data + "*.bam"):
        output_file_name = data_path + "/fastq_files/" + item.split("/")[-1].split(".bam")[0] + ".fastq"
        if os.path.exists(output_file_name):
            continue
        os.system("samtools fastq --threads 16 " + item + " > " + output_file_name)


def demultiplex_cutadapt(data_path):
    cut_command = "cutadapt -a file:barcoding_primers.fasta "
    flags = "--revcomp -m 1000 --action=lowercase -q 10 --cores=0 "
    output_path = "-o " + data_path + "{name}.fq "
    input = "all_reads.fastq"

    os.system(cut_command + flags + output_path + input)


def fastq_to_fasta(data_path):
    for file_name in glob.glob(data_path + "M*.fq"):
        os.system("sed -n \'1~4s/^@/>/p;2~4p\' " + file_name + " > " + file_name.split(".fq")[0] + ".fa")


def tiling_sequences(data_path, tile_len=20):
    tile_dict = {}
    counter = 0
    for fasta_file in glob.glob(data_path + "M*.fa"):
        with open(fasta_file, "r") as input_fasta:
            for line in input_fasta:
                if line.startswith(">"):
                    pass
                else:
                    tiles = [line.strip().upper()[0 + i:tile_len + i] for i in range(0, len(line.strip()), tile_len)]
                    counter += 1
                    for tile in tiles:
                        try:
                            tile_dict[tile] += 1
                        except KeyError:
                            tile_dict[tile] = 1
        output_name = fasta_file.split(".fa")[0] + "_unique_tiles.fa"
        with open(output_name, "w") as output_fasta:
            i = 1
            for entry in tile_dict.keys():
                i += 1
                output_fasta.write(">" + str(i) + "_" + str(tile_dict[entry]) + "\n")
                output_fasta.write(entry + "\n")


def hisat_index_build(data_path, fasta_file, index_name):
    os.system("hisat2-build " + fasta_file + " " + data_path + index_name)


def hisat_tile_alignment(data_path, index_name):
    for tile_file in glob.glob(data_path + "*_unique_tiles.fa"):
        header = tile_file.split("_unique_tiles.fa")[0]
        sam_output = header + "_output.sam"
        os.system("hisat2 -x " + data_path + index_name + " -f " +
                  tile_file + " -S " + sam_output + " --un " + header +
                  "_unaligned.fa --al " + header + "_aligned.fa --no-unal -p 16")


def find_tiles_unique_to_donorVSG(data_path):
    for sam_file in glob.glob(data_path + "*_output.sam"):

        align_dict = {}
        read_dict = {}
        with open(sam_file, "r") as input_sam:
            for line in input_sam:
                if line.startswith("@"):
                    continue
                else:
                    read_name = line.split("\t")[0]
                    read_dict[read_name] = line.split("\t")[9]
                    try:
                        align_dict[read_name].append(line.split("\t")[2])
                    except KeyError:
                        align_dict[read_name] = [line.split("\t")[2]]

        donor_aligned_reads = sam_file.split("_output")[0] + "_donor.fa"
        with open(donor_aligned_reads, "w") as output_fasta:
            for item in align_dict.keys():
                if "VSG4806" in align_dict[item] and "VSG8" not in align_dict[item]:
                    output_fasta.write(">" + item + "\n")
                    output_fasta.write(read_dict[item] + "\n")


def find_mosaic_reads(data_path, tile_align_threshold=50, tile_len=20):

    for donor_tiles in glob.glob(data_path + "*_donor.fa"):
        tiles = []
        with open(donor_tiles, "r") as input_tiles:
            for line in input_tiles:
                if line.startswith(">"):
                    pass
                else:
                    if len(line.strip()) == tile_len:
                        tiles.append(line.strip())

        file_name = donor_tiles.split("_donor")[0] + ".fa"
        output_name = donor_tiles.split("_donor")[0] + "_mosaic.fa"
        if os.path.exists(output_name):
            continue
        read_count = 1
        with open(file_name, "r") as reads, open(output_name, "w") as output:
            for line in reads:
                if line.startswith(">"):
                    read_name = line.strip()
                    print(file_name + str(read_count))
                    read_count += 1
                else:
                    count = 0
                    for entry in tiles:
                        if len(line.strip().split(entry)) == 2:
                            count += 1
                            if count >= tile_align_threshold:
                                output.write(read_name + "\n")
                                output.write(line)
                                break


def read_len_filter(data_path):
    for fasta_file in glob.glob(data_path + "M*.fa"):
        updated_name = fasta_file.split(".fa")[0] + "_2kb.fa"
        if "unique_tiles" in fasta_file or "aligned" in fasta_file or "donor" in fasta_file:
            continue
        with open(fasta_file, "r") as input_file, open(updated_name, "w") as output_file:
            for line in input_file:
                if line.startswith(">"):
                    read_name = line
                else:
                    if len(line.strip()) <= 2000:
                        output_file.write(read_name)
                        output_file.write(line)


def make_single_line(file_path):
    with open(file_path, "r") as fasta, open(file_path.split(".fa")[0] + "_ol.fasta", "w") as output:
        seq = ""
        name = ""
        for line in fasta:
            if line.startswith(">"):
                if seq != "":
                    output.write(name)
                    output.write(seq + "\n")
                    seq = ""
                name = line
            else:
                seq = seq + line.strip()


def seq_comparison(target, read, donor):
    """ Identify events within mosaic reads
    Target, donor, and the read are directly compared base by base. Random point mutations/errors are
    identified. Additionally, cross over events are identified as the range of identical bases between
    target and the donor where the read switches patterns from one to the other. Cross overs in any direction
    can be identified. A cross over event is only logged if the read takes advantage of two consecutive
    opportunities to switch. A list of all events found is returned.
    Args:
        target: string, the predicted sequence of target corresponding to the read
        read: string, the consensus sequence of the read, has some lower case information
        donor: string, the predicted sequence of the donor VSG corresponding to the read

    Returns:
        cross_over_record: a list of tuples, (start, end, event type) - all values correspond to locations
        within the read
            types of events:
            target_to_donor: cross over event from target to donor
            donor_to_target: cross over event from donor to target
            random_donor: a single bp matching the donor strand in a large target section
            random_target: a single bp matching target in a large donor section
            random_error: a bp matching neither target or the donor
            target_to_donor_edge: potential crossover event, but not enough proof of
                two consecutive matches with switch because of read length
            donor_to_target_edge: potential crossover event, but not enough proof of
                two consecutive matches with switch because of read length

    """

    # initialize variables
    f = "target"
    count = 0
    cross_over_record = []
    first = True
    previous_loc_donor = 0
    previous_loc_target = 0
    first_base_mismatch = False

    # loop through each sequence one bp at a time
    # read has some lowercase values in the primer
    for postarget, posread, posdonor in zip(target, read.upper(), donor):
        # store previous VSG
        previous_VSG = f

        # reset same to false
        # define the where the read matches first, the target or the donor
        same = False
        if count == 0:
            if posread == posdonor and posread != postarget:
                first_base_mismatch = True
                wrong_base = "donor"
                f = "donor"
                previous_VSG = "donor"
            elif posread != posdonor and posread == postarget:
                first_base_mismatch = True
                wrong_base = "target"
        elif first_base_mismatch:
            if wrong_base == "donor":
                if posread == posdonor and posread != postarget:
                    first_base_mismatch = False
            elif wrong_base == "target":
                if posread != posdonor and posread == postarget:
                    first_base_mismatch = False


        # identify type of event
        if posread != posdonor and postarget != posread:
            # read is unique, an error
            # store location of error in record
            random_event = (count, count, "random_error")
            cross_over_record.append(random_event)
        elif postarget == posread and posread != posdonor:
            # read matches target
            previous_loc_target = count
            f = "target"
        elif postarget != posread and posread == posdonor:
            # read matches donor
            previous_loc_donor = count
            f = "donor"
        else:
            # read, donor, and target are the same
            same = True

        # Ensures that cross over events take advantage of two consecutive opportunities to match
        # new fragment and not old fragment, bases that are all the same do not count
        if not first and not same:
            # if previous VSG fragment is the same as this one, this is a true event
            if f == previous_VSG and possible_range[0] != 1:
                cross_over_record.append(possible_range)
            elif first_base_mismatch:
                random_event = [(0, 0, "random_" + wrong_base)]
                cross_over_record = random_event + cross_over_record
                first_base_mismatch = False
                previous_VSG = f
            elif possible_range[0] <= 1 and possible_range[1] >= 1:
                previous_VSG = f
            else:
                # the previous VSG fragment is not the same, this is a single base cross over
                # These are counting as random and not events
                # position of the mismatch is reported
                random_event = (possible_range[1], possible_range[1], "random_" + previous_VSG)
                cross_over_record.append(random_event)
                # reset previous VSG as this is not a true cross over event
                previous_VSG = f
            # reset flag for next cross over
            first = True

        if f != previous_VSG:
            if f != "target":
                possible_range = (previous_loc_target + 1, count, previous_VSG + "_to_" + f)
            else:
                possible_range = (previous_loc_donor + 1, count, previous_VSG + "_to_" + f)
            first = False

        count += 1

        # end case; There are some reads where the final opportunity for a cross over event seems like
        # one, but there is no subsequent opportunity to check due to the length of the read (all bases are
        # the same).
        if count == len(read):
            try:
                if possible_range not in cross_over_record:
                    if len(cross_over_record) != 0:
                        if possible_range[0] > cross_over_record[-1][0]:
                            # any added range needs to be completely after any previously added records
                            end_info = possible_range[2] + "_edge"
                            end_cross_over = (possible_range[0], possible_range[1], end_info)
                            cross_over_record.append(end_cross_over)
            except NameError:
                pass
    return cross_over_record


def find_insertion_sites(data_path):
    with open(data_path + "all_corrected_mos_ol.fasta", "r") as input_fasta, \
            open(data_path + "all_mos_cross_overs.tsv", "w") as output_file:
        count = 0
        for line in input_fasta:
            if line.startswith(">"):
                read_name = line.strip()
            else:
                count += 1
                seq = line.strip()
                potential_20bp_list = ["GGAAAACCAAAACGATTCAA", "GGAAAACCAAAACGATTCCA", "GGAAAACACAAACGATACCA",
                                       "TGAAAACCAAAACGATTCCA",
                                       "GGAAAACCAAAACGTTCTCA", "GGAAAACCGAAACGATTCGA", "TTGTGACACACTTCATTAGT",
                                       "GGAAAACTGAGTCATTCCAA",
                                       "CCGAAACCGAAACGATTCCA", "CCGAAACCGAAACGATTCCA", "GGAAAACCGAAACTATTCCA",
                                       "GGAAAACCAAAACGATTCCA",
                                       "GGAAAACCGAAAATGTTGTG", "ACCAAACAATTCCAAACACT", "GGAAAACCGAAACCATTCCA",
                                       "GGAAAACCAAAACCATTCCA", "TCAAAACCGACTCGATTCCA", "GGAAAACCGAGACCATTTCG"]
                trimmed_seq = ""
                try:
                    trimmed_seq = "GGAAAACCGAAACGATTCCA" + seq.split("GGAAAACCGAAACGATTCCA")[1][:1552]
                except IndexError:
                    try:
                        trimmed_seq = seq.split("TACAGTTTCTGTACTATATT")[1][:1572]
                    except IndexError:
                        try:
                            trimmed_seq = "GGAAAACCGAAACGATTCCA" + "AACACTATTCGACGTTAAGC" + \
                                          seq.split("AACACTATTCGACGTTAAGC")[1][:1532]
                        except IndexError:
                            for bp20 in potential_20bp_list:
                                try:
                                    trimmed_seq = bp20 + seq.split(bp20)[1][:1552]
                                    if len(trimmed_seq) == 1572:
                                        break
                                except IndexError:
                                    pass

                if len(trimmed_seq) != 1572:
                    print(read_name)
                    pass

                record = seq_comparison(vsg_8_seq, trimmed_seq, vsg_4806)
                if len(record) >= 300:
                    pass

                for item in record:
                    output_file.write(
                        "\t".join([read_name.split()[0].split(">")[1], str(item[0]), str(item[1]), item[2]]) + "\n")


def read_name_leaders(data_path):
    first = True
    # get consol_groups
    with open(data_path + "all_mos_clusters_99.fa.clstr", "r") as input_cluster_file:
        clusters = {}
        read_list = []
        for line in input_cluster_file:
            if "Cluster" in line:
                if not first:
                    clusters[leader_name] = read_list
                read_list = []
                first = False
            elif "at" not in line:
                leader_name = line.split(">")[1].split("...")[0]
                read_list.append(line.split(">")[1].split("...")[0])
            else:
                read_list.append(line.split(">")[1].split("...")[0])
        clusters[leader_name] = read_list

    invert_clusters = {}
    for entry in clusters:
        for item in clusters[entry]:
            invert_clusters[item] = entry

    sample_info_dict = {}
    for file_name in glob.glob(working_path + "M*_mosaic.fa"):
        sample_info = file_name.split("_mosaic")[0].split("/")[-1]

        read_count = 1
        with open(file_name, "r") as reads:
            for line in reads:
                if line.startswith(">"):
                    read_name = line.strip().split()[0]
                    read_count += 1
                else:
                    sample_info_dict[read_name] = sample_info

    with open(data_path + "read_name_sampleinfo.tsv", "w") as output_tsv:
        for item in invert_clusters:
            print(item)
            output_tsv.write("\t".join([item, sample_info_dict[">" + item], invert_clusters[item]]) + "\n")


def main():
    work_path = ""
    convert_bam_to_fastq(work_path)
    os.system("cat " + work_path + "fastq_files/*.fastq > " + work_path + "all_reads.fastq")
    demultiplex_cutadapt(work_path)
    fastq_to_fasta(work_path)
    tiling_sequences(work_path)
    hisat_index_build(work_path, "vsg8_donor.fasta", "vsg8_donor")
    hisat_tile_alignment(work_path, "vsg8_donor")
    find_tiles_unique_to_donorVSG(work_path)
    find_mosaic_reads(work_path)
    os.system("cat " + work_path + "*_mosaic.fa > " + work_path + "all_mosaic.fa")
    read_len_filter(work_path)
    os.system("cd-hit-est -i " + work_path + "all_mosaic.fa" + "-o" + work_path + "all_mos_clusters_99.fa -c 0.99 -n 8 -M 16000 -d 0 -T 16")
    make_single_line(work_path + "all_mos_clusters_99.fa")
    # these were then corrected by hand into all_corrected_mos_ol.fa
    find_insertion_sites(work_path)
    read_name_leaders(work_path)


if __name__ == '__main__':
    main()

