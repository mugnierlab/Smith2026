"""
From Beaver et al. find the candidate mosaic recombination events
Smith 2024
Author: Jaclyn Smith
"""

antat_seq = "ATGGTCACCAAGGAGCGAAACGCAGCATTAAAAATTGTAATGTTAGTCGCTTCAGCACTGACACTACACCCACAACAAGCTCTAGCTCAGACCGCTGGTAGGCCCCTTGCAGATGTGGTAGCCAAAACTCTATGTACTTATTCAAAAACGGCCAAACGCCAGGCAGCAAACCTGGCGCAAACACTACAACGAGCCAGCTCAGCAGCAAAGCAATCCAGACAAGCGCAGCAGTTAGCGGCTTTAGCACTGGCCAAACTGCCAGACTACAAAGAAGCAGCCGCGACACTGTTAATTTACGCCACGCACAAAATACAAGACGCGCAAGCCAGCATCGAAAACTGGACAGGAGAGAATACTAAGCTAGTTGGCCAGGCGATGTATTCCTCAGGGAGAATCGACGAACTGATGTTGCTACTAGAAGGGCACCGAGAGGACGGCGCGAACGGACAGGACAAAACTTGCCTAGGCGCGGCCGCCGGCGGCAATACAGTAAATGAATTCGTCAAAACAGAATGCGACACGGAAAGCGGCCACAACATCGAGGCAGACAACTCAAACATAGGGCAAGCGGCAACGACTCTAAGCCAAGAAAGTACAGACCCAGAAGCCAGCGGAGGCGCAAGCTGCAAAATAACAGCAAACCTTGCCACTGACTACGACAGCCATGCGAATGAGTTACCGCTACTCGGCGGCCTGCTAACCATACACAACGCAGGCGGCTTCAAAACAGGACAAAGCTTGCAAACCGCAGCACCAACCAACAAGCTAATCAGCGCACTCAAAAATAAGGGCGCCGGTGTCGCAGCTAAACTGGCAACTGTAACGTCGGCAGCACCTACAAGCAAGCAGGAACTCAAAACACTACTGGCTTCGAAAGGGGAACGCGCCAAACTCCAAGCAGCCAACGACGAGTATAATAACTGGAAACCAGGCGCCAAGCCTGAGGACTTCGACGCCCACATCAAGAAAGTGTTCGGCGCAGAAGACGGCAAAGACAGCGCCTATGCCATTGCACTTGAAGGAATATCCATTGAGGTTCCCCTCGGAGGAGGACAAACACAAAACAAACAACTCTATTCCATGCAGCCAAAAGACCTAATGGCAGCTTTAATAGGAACGATAGCAGAACTCCAAACAGCCGCAGCAACCAAACCAGCATGCCCAGGCCATAAACAAACAACCACGGAAAGTGACGCCCTATGCAGTAAAATAAAGGATGCAAACGAATGCAACAGCAAGCCTTACTGCAGTTATAACGAAACCGCAGCTTATGGCGACAAAAAGTGCCAATTTAATGAAACAAAGGCCTCAAAAAATGGAGTCCCTGTAACGCAAACTCAAACTGTAGAACCTTCAACCAATCCAGAAAAGTGCAAAGGGAAAGAAGAGAAAGACTGCAAATCCCCGAATTGTAAATGGGAGGGCGAAACTTGCAAAGATTCCTCTATTCTAGTAAACAAGAAATTCGCCCTCAGCGTGGTTTCTGCCGCATTTGTGGCCTTGCTTTTTTAATTTTCCCCCTCTTTTTCTTGCTAAAAATTCTTGCTACTTGAAAAATTTTCTGATATATTTTAACAC"
donors = {
    "Tbb1125VSG-228": "ATGGTCGCTAAGAAGTGCAGCGCAGCATTAAAGATAGTAATGTTAGTCGGTGCCGCACTGACACTACACCAACAACAAGCTCTAGCTCAGACCGCTGGTAGGCCCCTTGCAGATGCGGTAGGCAAAGCACTCTGCACTTATTCAAAGACAGCCAAACGACAGGCAGCAAACCTAGCCCAAGCTCTAGATCGCGGCATCACAGCAGCAAAAAAGTCGCAACAAGCGCAGCAGTTAGCGACGATAGCACTGGCGAAACTACCCCACTACAGAGAAGCAGCAGCGACGATCCTCATTTACGCCAAAAACAAAAGAGCAGAAGCAGAAGCAAACATCGAAAACTGGAAAGGCCAGAAAACCAAACTGGTGGGGCAGGCAATGTATTCCTCAGGCAGAATCGACGAGCTGATGTTGATGCTAGAAGGCCACAGAGACGGACAATCAGCAGGGCAAACCAAAACTTGCCTAGGCGCAGCAGGAAACGGCAACACAGTAGATGAATTCGTCAAAACAGAGTGCGACACGGAACAAGACCACAACATCAACGCCGACGAATCAGACATAGAACAAGCGGCAGCAACCCTAAGCCAGGAAAATAGAGACCCGGAAGCAGGCGGCGGAACCAACTGCAAGATCACAGGCAACCTTGCCAGCGACTACGACAGTCATCCAAATGATCTGAGTTTGCTAGGCGGACTGCTAACAATACACAATGGGGGCGGCTTTAAGGCAACGACAACAATAAAAACCGCGGCGGCGGGCAACAAACTAATCAGCGCCCTCGCGAGCAAGGTTAACGACATCGCTGCTAACCTCAAAGCACACACGGAATCGGCACCAACGACCAAACAAGAACTCAAGACACTACTGGGCAGCAAAGGAGCACGGAGCAAGCTAGAAGCAGCCAACGACGAGTACAATAGCTGGGAAGCAGGAAAAAAGCCTGTAAACTTCGACGAGCACATCAAAAAAGTGTTCGGCGCGGAAGACGGCAAAGACAGCGCCTATGCCCTTGCACTTGAAGGGATATCCATTGAAGTTCCACAAAAACCAGGGACCACAGAAAGCAAACAACTCTATTCCATGCAGCCAAAAGACCTAATGGCAGCTTTAATAGGAACCATAGCAGAAATACAAAAAGCCGCAGCGACAAAAGCACCATGCCCAAAGCATAAACTGACAAGCGCTGAAAGTGACGCCCTATGCAGTAAAATAAAGGATGCAAACGAATGCAACAGCAAGCCTTTCTGCAGTTATAACAGCACCGAAACTGACACAGCTAAAAAGTGCCAATTTAATGAAACCAAAGCTGACAAAAGTGGAGTTTCGTTGCCTAAAACCGGACCTACCGGTACTGAAGCAACAACTGATAAGTGTAAAGATAAGACAAAAGATGAGTGCAAATCTCCGAATTGTAAATGGGAGGGCGAAACTTGCAAAGATTCCTCTATTCTAGTAAACAAACAATTCGCCCTCAGCATGGTTTCTGCTGCCTTTATCGGTTTGATATCATTTTAGAATTGTAAGGATTTTTATAAATTTTATGAGTTACGAAAACTTATTATTTTGAGAAACTTCTTAATACTTGATATAATTTAATATATTTTAACATTTTTCGGTAGAAATGTGCGAATTAGCGATTATACAAATGTTGGAAAATATAAATAAAATG",
    "Tbb1125VSG-3110": "ACAAGCTCTAGCTCAGACCGCTGGTAGGCCCCTTGCAGATGTGGTAGGCAAAACTCTATGTACTTATTCAAAAACGGCCAAACGCCAGGCAGCAAACCTGGCACAAACACTAGAACGAGCCAGCTCAGCAGCAAAGCAATCCAGACAAGCGCAGCAGTTAGCGGCTTTAGCCCTGGCGAAACTCCCAGAGTACAAAGAAGCAGCCGCGACGCTGTTAATTTACGTCATGATGAAAGTGGAAGCAGCACAAGCAAGCATCGAAAACTGGACAGGAGAGAAGACAAACTAGTAGAGCAGGCGATGTATTCCTCAGGCAGGATAGACGAGCTGATGTTGCTCCTAGAAGGGCACCGGGAGGACGGACCGAGCGGTGCGGACAAAACTTTCCTAGGCGCGGCCGCCGGCGGCAATACAGTAAATGAATTCGTCAAAACAGAATGCGACACGGAACACGACCACAACATCGAGGCAGACAACTCAAACATAGGGCAAGCGGCGACGACTCTAAGCCAAGAAAGTACAGACCCAGAGGCAGCCGGCAGAAGCGACTGCAAGGTCACAACCAACCTTGCTAGCGACTACGACAGCCATGCGAATGAACTACCGCTACTCGGCGGACTGCTAACCATACACAATGGGGGCGGCTTCAAAAGCGGTCAGACCATACAAACTACAGCAACCACGAACAAGCTAATCAGCGCACTCAAAAATAAGGGCGCCGGTGTCGCCGAGAGCCTAAAAACTGTAACAGCGGCAGCACCTACAAACAAACAGGAGTTCAAGACACTACTGGCCTCAAAAGCCGAGCGTGCCAAACTGCAAGCAGCGAACGAAGAGTATAATAACTGGAAACCAGGCGCCAAGCCCGCGGACTTCGATGCCCACATCAAAAAAGTGTTCGGCGCGGAAGACAACAAAGACAGCGCCTATGCCCTTGCACTTGAAGGGATATCCATTGAGGTTCCCCTCGGAGGAGGACAAACACAAAACAAACAACTCTATTCCATGCAGCCAAAAGACCTAATGGCAGCTTTAATAGGAACGATAGCAGAACTCCAAACAGCCGCAGCGACGAAACCAGCATGCCCAGGCCATAAACAAACAACCACTGAAAGTGACGCCCTATGCAGTAAAATAAAGGATGCAAACGAATGCAACAGCAAGCATTTCTGCAGTTATAACGGCACCGAAACTGACTCAGCTAAAAAGTGCAAATTTAATGCCTCAAAAGCTTCAGCAAGTGGTGCCCCTGTAACACAAAGTCAAACCGGAGGGAGTGAAACAACAACAGAGAAATGCAAAGGCAAGAAAAAGGATGACTGCAAATCTCCGGATTGTAAATGGGAAGGCGAAACTTGCAGAGACTCTAGTTTCTCATAAATAAGAAATTTTCTCTGATAGCTGCTGCTTTGTGAATTTGGTGGTATTTTAGAATTTTAATGATTTTTATGCAATTTTATGAAAATTTATGAAAATTTTCTATTTCGAGATAATTTAAAATATTTAAGCACTTTATGTTGGAAATGTGAAAATGAATAAAAATATGGAAATACAAAAAATATGAAATTTTGGTATGATAGTGAACTATGAGTGCAAAATAGTAGGACTACAGGATGTGAAAAAATAGCGCTATTACATAATTTGGTTTATTTTTTTGTTTGATGAATATTAAAAGCGCTT",
    "Tbb1125VSG-2986": "ATGGTCGCCAAGGAGTGCAACGTAGCATTAAAAATCGTAATGCTAGTCGCTTCAGCACTGACACTACACACACAACAAGCTCTAGCTCAGACAGCTGGTAGGCCCCTTGCAGATGTGGTCGGCAAAGCTCTATGTACTTATTCAAAAACGGCGAAACGACAGCCAGCGAATCTAGCCCAAGCTCTAGATCGCGCCATCACAGCAGCAAAGAAATCGGAGCAAGCACAGGCGTTAGCAGCCGTTGCCCTGGCCAAGCTCCCAGACTATCAGCAAGAAGCCGGAACGCTGTTAATTTACGTCAGGATGAAAGTGGAAGCAGCACAAGCGAGCATCGAAAACTGGACAGGAGAAAAGACAAAACTGGTAGAGCAGGCGATGTATTCCTCAGGCGGGATAGACGAGCTGATGTTGCTCCTAGAAGGGCACCGGGAGGACCGGGCGAACGGCGCCGACAAAACCTGCCTAGGCGCAGCAACAAGCGGAAACACAGTCGCTGAATTCGTCAAAAGAAATTGCGATACGGAAAATGATCACGACATCAGCGCGGACAACTCAGACATAGGCCAAGCGGCGGCAGATCTAAGCCAAGCAAGCACAGACACAGAGGCAGGCGGCGGAAGCAACTGCAAAATAACAGCAAACCTTGCCAGCGACTACGACAGTCATGCGAATGAACTAACGCTACTCGGCGGGCTGCTAATCATACACAATGAAGGTGGCTTTAAAAGCGGTCAAACCATAAAAACCGCAGCAACCGCGAACAAGCTAATCAGCGCCCTCAAAAATAAAGGAACCAGCGTTGCCAAGAGCCTAAAAACAGTAACAGCGGCAGCACCTACAAACAAACAGGAGTTCAAGACACTACTGGCTTCGAAAGGGGAGCGTGCAAAACTCCAAGCAGCCAACGACGAGTATAATAACTGGAAACCAGGCGCCAAGCCTGTGGACTTCGATGCCCACATCAAGAAAGTGTTCGGCGCAGAAGACGGAAAAGACAGCGCCTATGCCATTGCACTTGAAGGGATATCCATTGAGGTTCCCCACGGAGCAAGGAAAACAGAAAGCAAACAACTTTATTCCATGCAGCCAAAAGGCCTAATGGCAGCTTAAATAGGAACCATAGCAGAACTTCAAACAGCCGCAGCGACGAACACAGTTTGGCCAAAACATAAACAAACAACCACTGAAAGTGACGCTTTATGCAGTAAAATTAAGGACGCAAACAAATGCAGCAGCAAGCCTTTCTGCAGTTATAACGGCACAGAAACTGAATCAGCTAAAAAATGCAAATTTAATGCGACAAAAGCTGAAAACTCCGGCGCTTCTGTAACACAAAGTCAAACCGGAGGGAGTGAAACAACAACAGAGAAATGCAAAGGCAAGAAAAAGCAGGAGGACTGCAAATCTCCG",
    "Tbb1125VSG-7358": "GCGGCTCTACCGCTACTCGGCGGACTGCTAACCATACACAATGGGGGCGGCTTCAAAAACGGCCAGACCATACAAACCACAGCAACCACGAACAAGCTAATCAGCGCACTTAAAAATAAAGGAGCCGGTGTCGCCGAGAGCCTAAAAACTGTAACATCGGCAGCACCTACAAACAAGCAGGAGCTCAAGACACTACTGGCTTCAAAAGGGTAGCGTACCAAACTCCAAGCAGCCAACGACGAGTATAATAACTGGAAACCAGGCGCCAACCTTGAGGACTTCGATGCCCACATCAAGAAAGTGTTCGGCGCAGAAGACAGCAAAGACAGCGCCTATGCCATTGCACGTGAAGGGATATACATTAAAGTTCCAACTAAAGCAGGGAAAACAGAAAGCAGACAACTCTATTCCATGCAGCCAAAAGACTTAATGACAGCTTTAATAGGAACCATAGCAGAACTCCAAACAGCCGCAGCGACGAAACCAACATGCCCAGACCATAAACGGGCAACCACTGAAAGTGACGCCCAATGCAGTAAAATAAAGGATGAAAACGAATGCAACAGCAAGCATTTCTGCAGTTATAACGGCACCGAAACTGACTCAACTAAAAAATGCAAATTTGATGCCAAAAAAGCCACATCAAATAGTGTCTCTGTAACACAAACTCAAACTGGAGGACAGAAACTGCAACAGGAAATCGCAAAGGGAAAGAACAAAAAGACTGCAAATCTGCGGATTGTAAATGGGAGAGAACAGATTGCAAAGATTTCAGCATTCTCGTCAATAAGAAATTGGCTCTGAGTATGGCTGCTGCTTTTGTTAGTTTTGCGGCAATTTAAAATTCTATAATTTTAAGGATATTTGCTCAATTTTGCAATATTTATGAAATTTCTTGTTTTGGGAGAATTTGCTAACATATGATCTATTTTAACAACTTTTTGTCGGATATGTGAAAATAAGTAAAAATATGGGAATACGAAAAATATGAAATTTTGGTATAATAGAGGACTATTGGTGCAAAATAGTGGGGGTCAACGTATGGAAAATCAACTATAATATAATTTGGCTTATTTTGTTTGTTTGATGAATATTAAAAGCGCT"
}

donor_starts = {
    "Tbb1125VSG-228": 0,
    "Tbb1125VSG-3110": 74,
    "Tbb1125VSG-2986": 0,
    "Tbb1125VSG-7358": 676
}


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
            # elif possible_range[0] <= 1 and not first_base_mismatch:
            #     # in the case of read1, there are some instances where the read is entirely donor, but parts
            #     # match target - don't want to confuse those parts for a cross over.
            #     # THis makes sense, but I will remove this logic and put it to test whether we should use this function
            #     #come back and take this off once completed !!!!!!@@#@#$@$@#@$$
            #     # single bp or no bp differences unique to target should not count as being a cross over event
            #     # now I will allow the random base to be executed here I think
            #     random_event = (possible_range[0], possible_range[1], "random_" + previous_VSG)
            #     previous_VSG = f
            # # elif possible_range[0] <= 1 and first_base_mismatch:
            #     random_event = [(0, 0, "random_" + wrong_base)]
            #     cross_over_record = random_event + cross_over_record
            #     previous_VSG = f

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
                            # if possible_range[0] == cross_over_record[-1][0] and \
                            #     possible_range[1] == cross_over_record[-1][1] and \
                            #     possible_range[2].split("_")[2] == cross_over_record[2].split("_")[1]:
                            # any added range needs to be completely after any previously added records
                            end_info = possible_range[2] + "_edge"
                            end_cross_over = (possible_range[0], possible_range[1], end_info)
                            cross_over_record.append(end_cross_over)
            except NameError:
                pass
    return cross_over_record


with open("check_for_mosaics.fa", "r") as names, open("all_orfs.fa", "r") as vsg_seqs:
    VSGs = {}
    for line in names:
        VSGs[line.strip()] = []
    first = True
    seq = False
    for line in vsg_seqs:
        if line.startswith(">"):
            if not first and seq:
                VSGs[name.strip()].append(sequence)
                seq = False
            name = line.strip()
            if name in VSGs.keys():
                first = False
                seq = True
                sequence = ""
        elif seq:
            sequence = sequence + line.strip()
    if not first and seq:
        VSGs[name.strip()].append(sequence)
count = 0
for item in VSGs:
    antat_start = "ATGGTCACCAAGGAGCGAAA"
    test_split = VSGs[item][0].split(antat_start)
    v228_start = "ATGGTCGCTAAGAAGTGCAG"
    v228_split = VSGs[item][0].split(v228_start)

    # these sequences start with antat!
    if len(test_split[0]) == 0:
        count += 1
        pass
        for donor in donors:
            start_pos = donor_starts[donor]
            VSGs[item].append((donor, seq_comparison(antat_seq[start_pos:], VSGs[item][0][start_pos:], donors[donor])))
    elif len(v228_split[0]) == 0:
        count += 1
        for donor in donors:
            start_pos = donor_starts[donor]
            VSGs[item].append((donor, seq_comparison(donors[donor], VSGs[item][0][start_pos:], antat_seq[start_pos:])))


with open("all_orfs.fa.clstr", "r") as cluster, open("multi_crossOvers_all_WTmice_unedited.txt", "r") as mosaic_records:
    master_list = {}
    leader_list = {}
    list_to_append = {}
    first = True
    for line in mosaic_records:
        info = line.replace("'", "").split(".")
        master_list[info[0]] = line.replace("'", "")
    for line in cluster:
        if "Cluster" in line:
            if first:
                group = []
                leader = ""
                first = False
            else:
                leader_list[leader] = group
                # reset stats
                group = []
                leader = ""
        elif "*" in line:
            leader = line.strip().split("...")[0].split(", ")[1]
        else:
            group.append(line.strip().split("...")[0].split(", ")[1])
    leader_list[leader] = group
    for item in master_list.keys():
        for entry in leader_list[item]:
            list_to_append[entry] = master_list[item]
with open("abex7_all_orfs.fa", "r") as vsg_seqs:
    VSGs_check = {}
    for line in list_to_append:
        VSGs_check[line] = []
    first = True
    seq = False
    for line in vsg_seqs:
        if line.startswith(">"):
            if not first and seq:
                VSGs_check[name.strip()].append(sequence)
                seq = False
            name = line.strip()
            if name in VSGs_check.keys():
                first = False
                seq = True
                sequence = ""
        elif seq:
            sequence = sequence + line.strip()
    if not first and seq:
        VSGs_check[name.strip()].append(sequence)
    with open("check_leng.fa", "w") as output2:
        for item in VSGs_check:
            output2.write("%s\n%s\n" % (item, VSGs_check[item][0]))
    with open("check_leng.fa", "a") as output2:
        for item in VSGs:
            output2.write("%s\n%s\n" % (item, VSGs[item][0]))
    list_to_append['>TRINITY_DN79_c0_g1_i2_1_1_911_RC_M30_D10_Ear'] = '>TRINITY_DN6_c0_g1_i5_1_0_1053_RC_C3M4_D14_subcu_R1.(895, 914, donor_to_target, VSG-228)\n'
    list_to_append['>TRINITY_DN29_c0_g1_i2_1_2_1420_RC_M35_D14_SubCu_R1'] = '>TRINITY_DN76_c0_g1_i1_1_57_1528_M3_D14_blood_R1.(1185, 1243, target_to_donor, VSG-228)\n'
    master_list.update(list_to_append)
    master_list.update(list_to_append)

with open("single_COs.csv", "w") as output:
    count = 0
    for item in master_list:
        info = master_list[item].split(".")
        try:
            sample_info = info[0].split("_M")[1].split("_")
            mouse = "M" + sample_info[0]
        except IndexError:
            sample_info = info[0].split("_C")[1].split("_")
            mouse = "C" + sample_info[0]
        print(item)
        day = sample_info[1]
        genotype = sample_info[2]
        # placeholders
        primer = "1R"
        direction = "R"
        for entry in info:
            if ">" in entry:
                pass
            else:
                count += 1
                co_info = entry.split(", ")
                output.write(", ".join([mouse, day, genotype, primer, direction, str(count), co_info[3].split(")")[0],
                                        str(86), co_info[0].split("(")[1], co_info[1], co_info[2], "\n"]))


    print(line)


print("test")
