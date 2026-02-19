"""VSG-AMP-Seq: This module identifies true mosaic events from candidates. Recombination sites are
    identified
   Smith 2024
   Author: Jaclyn Smith
"""
import Levenshtein as lev
import regex
import global_target as gbl
import aux_functions


def target_primer_id(primer_dict):
    """
    determines ranges of potential read2s from primer sequences
    Args:
        primer_dict: dictionary from sort_reads_by_primer key = sequence, value = primer_name
    Returns:
        primer_ranges: a dictionary, key = primer name; value = tuple with primer range
    """
    primer_ranges = {}
    allowed_mismatch = str(0)
    for entry in primer_dict:
        primer_name = primer_dict[entry]
        primer_seq = entry
        primer_direction = primer_name[1]
        if primer_direction == "F":
            primer_regex = "(?b)(" + primer_seq + "){e<=" + allowed_mismatch + "}"
            primer_return = regex.search(primer_regex, gbl.full_length_target)
            start = primer_return.span()[0]
            end = start + 150
        else:
            cor_seq = aux_functions.rev_complement(primer_seq)
            primer_regex = "(?b)(" + cor_seq + "){e<=" + allowed_mismatch + "}"
            primer_return = regex.search(primer_regex, gbl.full_length_target)
            end = primer_return.span()[1]
            start = end - 150

        primer_ranges[primer_name] = (start, end)
    return primer_ranges


def read_in_donorVSGs(vsg_file):
    """
    This function reads in the VSG donor FASTA file into a dictionary
    Args:
        vsg_file: the file name ending in .fasta

    Returns:
        vsg_candidate_dict: dictionary - key: VSG_name, value: VSG sequence

    """
    # look for reads aligning to alternative VSGs with primer sequence
    vsg_candidate_dict = {}
    with open(vsg_file, "r") as VSGs:
        title = ""
        sequence = ""
        for line in VSGs:
            if ">" in line:
                title = line.strip().split(">")[1].split(" : ")[0].split("Tbb1125")[1]
            else:
                sequence = line
            if title != "" and sequence != "":
                vsg_candidate_dict[title] = sequence

                title = ""
                sequence = ""
    return vsg_candidate_dict

# just keep in mind the gen code here was different as far as making everything uppercase


def rev_complement(read):
    """
    create reverse complement of DNA sequence
    Args:
        read: string, DNA sequence to be converted

    Returns:
        rev_comp_read: string, converted DNA sequence

    """
    comp_read = ""
    for position in read:
        # gen_code is a dictionary, lowercase letters are replaced with uppercase versions
        new_base = gbl.gen_code[position]
        comp_read = comp_read + new_base
    # reverse the read and return
    rev_comp_read = comp_read[::-1]
    return rev_comp_read


def get_donor_range(consensus, un_anch_seq, VSG_name, donor_VSGs):
    """ Range for donor VSG overlapping the read
    This uses the short seq identifying donor VSGs to locate the region with the donor matching the full
    length of the read.
    Args:
        consensus: string, the consensus sequence from a combination of Read1 and Read2
        un_anch_seq: string, the short donor VSG sequence used to identify the donor VSG(s), info[1]
        VSG_name: string, the name of the donor VSG to find the range with

    Returns:
        donor_start: int, start position in donor VSG that aligns with read
        donor_end: int, end position in donor VSG that aligns with read
        indel: string, a string indicating types of errors found in un_anch_seq compared to the donor
            This will not match the needed length and will need to be taken into account later
    """
    # check that unanchored sequence matches consensus, with quality score dictating sequence, sometimes un_anch_seq no
    # longer found
    if consensus.count(un_anch_seq) != 1:
        # allowing 2 mismatches at the moment
        allowed_mismatch = str(2)
        print("un_anch_seq doesn't match consensus, or consensus messed up")
        un_anch_regex = "(?b)(" + un_anch_seq + "){e<=" + allowed_mismatch + "}"
        un_anch_loc = regex.search(un_anch_regex, consensus)

        if un_anch_loc:
            un_anch_seq = consensus[un_anch_loc.span()[0]: un_anch_loc.span()[1]]
        else:
            print("This consensus read significantly different from R1 sequence."
                  "Consensus needs to go back through mosaic discovery process")
            return 0, 0, ""

    # divides consensus into two parts, separated by previously identified donor VSG piece
    read_parts = consensus.split(un_anch_seq)
    # initialize variables
    allowed_mismatch = str(1)

    # donor VSG regex, this returns a search object that identifies the location within the donor where
    # the match occurs
    donor_un_anch_regex = "(?b)(" + un_anch_seq + "){e<=" + allowed_mismatch + "}"
    donor_vsg_loc = regex.search(donor_un_anch_regex, donor_VSGs[VSG_name])
    if donor_vsg_loc:
        # to handle instances where the location of un_anch_seq ends up with a deletion
        # and the found insert is too short
        if len(un_anch_seq) < donor_vsg_loc.span()[1] - donor_vsg_loc.span()[0]:
            indel = "insertion"
        elif len(un_anch_seq) > donor_vsg_loc.span()[1] - donor_vsg_loc.span()[0]:
            indel = "deletion"
        else:
            indel = "none"

        # using the position of the donor VSG piece and the consensus sequence,
        # the start and end positions are calculated
        donor_start = donor_vsg_loc.span()[0] - len(read_parts[0])
        if len(read_parts) != 1:
            donor_end = donor_vsg_loc.span()[1] + len(read_parts[1])
        else:
            donor_end = donor_vsg_loc.span()[1]
        return donor_start, donor_end, indel
    else:
        return 0, 0, ""


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


def narrow_down_donor(donor_co_record):
    """
    Removes donors where a better match is clearly present. First, I detect if there is a mismatch
    present in the cross over record. If absent from at least one record, those without mismatches
    are identified and returned. otherwise no record changes are made.
    Args:
        donor_co_record: a list, entries into the list are organized by donor VSG, identified events
        and the target start

    Returns:
        likely_donor: a list, entries into the list are organized by donor VSG, identified events
        and the target start. trimmed so just the more likely donors are returned
        donor_co_record: a list, the input to the function

    """
    # check for randoms mismatches present in the alignments
    randoms_in_events = []
    for i, item in enumerate(donor_co_record):
        for event in item[1]:
            if "random" in event[2]:
                randoms_in_events.append(i)
                break
    # if a mismatch is not present in one of the alignments, that is the best alignment
    # any alignments without a mismatch are then returned
    if len(randoms_in_events) != len(donor_co_record) and randoms_in_events != []:
        likely_donor = donor_co_record
        for removal in reversed(randoms_in_events):
            likely_donor.pop(removal)
        return likely_donor
    # otherwise do not edit down the number of donors and return the list unedited
    else:
        return donor_co_record


def translate_protein(sequence, start, ATG_pos):
    """
    Protein translator given a DNA sequence
    This particular translator takes the start position, relative to the typical start of
    the protein in question.
    The reading frame is determined and sequence is looped through by 3s to obtain protein seq.
    Args:
        sequence: string, sequence to be translated
        start: int, position where substring starts in target sequence, to determine reading frame
        ATG_pos: int, position where coding sequence starts in target sequence

    Returns: string, translated protein sequence

    """
    # initialize variables
    protein = ""

    # determine reading frame
    # frames are 0, 1, 2
    rf = (start - ATG_pos) % 3

    # loop through sequence by 3, include final set with length + 1
    # rf 0 is included since it backtracks by 3
    for triple in range((3-rf), len(sequence)+1, 3):
        # only translate if length of triple is 3
        if len(sequence[triple-3:triple]) == 3:
            # some unknown codes (that have N) are translated as U for unknown
            # protein code key: triplet code
            # protein code value: single letter amino acid
            try:
                protein += gbl.protein_code[sequence[triple-3:triple]]
            except KeyError:
                protein += "U"
    return protein


def difference_check(string_1, string_2):
    """ Differences between two strings
    Determines number and position of identity differences between two strings
    This is generic and can work for proteins or DNA/RNA
    Args:
        string_1: string
        string_2: string

    Returns:
        count: int, number of differences
        changes: list of tuples, contains all the changes between the two strings

    """
    # using regex to search for matches between the two strings, only allowing mismatch errors
    difference_regex = "(?b)(" + string_1 + "){s}"
    difference_matches = regex.search(difference_regex, string_2)

    if difference_matches:
        changes = difference_matches.fuzzy_changes[0]
        count = difference_matches.fuzzy_counts[0]
    else:
        # major error!
        changes = []
        count = 0

    return count, changes


def identify_double_COs(cross_overs, target_pred_seq, read, start_pos, donor_VSG):
    """
    Identifying all double cross overs and related protein/DNA changes
    This function identifies start and end of a cross over event contained within a read
    The predicted target sequence and actual read sequence are translated and compared, the resutls
    are returned.
    Args:
        cross_overs:list of cross over events
        target_pred_seq: string, predicted sequence in target corresponding to read
        read: string, sequence from the read
        start_pos: int, start position of read within target
        donor_VSG: string, name of donor VSG

    Returns:
        individual_cross_overs: list of double cross over events and details

    """
    # start_pos is start or primer_range[0]
    double_cross_overs = []

    # sort out only events from the list
    pairs = list(filter(lambda events: events[2] in ["target_to_donor", "donor_to_target", "donor_to_target_edge"],
                        cross_overs))

    # for every two events, a double cross over,
    for pair in range(0, len(pairs), 2):
        if pair + 1 != len(pairs):
            pinsert_target = translate_protein(target_pred_seq[pairs[pair][0]: pairs[pair + 1][1]],
                                              start_pos + pairs[pair][0], gbl.protein_start)
            pinsert_read = translate_protein(read[pairs[pair][0]: pairs[pair + 1][1]],
                                                  start_pos + pairs[pair][0], gbl.protein_start)
            prot_ins_diff_count, prot_differences = difference_check(pinsert_read, pinsert_target)
            nt_ins_diff_count, nt_differences = difference_check(target_pred_seq[pairs[pair][0]: pairs[pair + 1][1]],
                                                                 read[pairs[pair][0]: pairs[pair + 1][1]])
            nt_length = len(target_pred_seq[pairs[pair][0]: pairs[pair + 1][1]])
            p_length = len(pinsert_target)
            if p_length > 0:
                p_ratio = prot_ins_diff_count / p_length * 100
            else:
                p_ratio = "undefined"
            nt_ratio = nt_ins_diff_count / nt_length * 100
            double_cross_overs.append((start_pos, pairs[pair], pairs[pair + 1], nt_ratio, nt_length, p_ratio, p_length,
                                       donor_VSG))
    return double_cross_overs


def primer_gap_test(consensus, primer):
    """
    Tests for a gap between the primer and the sequence, some of these had weird errors near primer.
    Args:
        consensus: a string, the consensus sequence in 5' to 3'
        primer: a string, primer name

    Returns:
        error string: a string, reports any error detected
        new_coord: an int, returns the new coordinates if error detected

    """
    if "F" in primer:
        primer_seq = consensus[21:]
    else:
        primer_seq = consensus[-21:]
    allowed_mismatch = str(2)

    target_regex = "(?b)(" + primer_seq + "){e<=" + allowed_mismatch + "}"
    target_match = regex.search(target_regex, gbl.full_length_target)
    if target_match:
        if target_match.fuzzy_counts[1] > 0:
            # this is an insertion by fuzzy changes, but is a deletion when compared
            # to the sequence
            if "F" in primer:
                new_coord = target_match.span()[1]
            else:
                new_coord = target_match.span()[0]
            return "deletion", new_coord
        elif target_match.fuzzy_counts[2] > 0:
            if "F" in primer:
                new_coord = target_match.span()[1]
            else:
                new_coord = target_match.span()[0]
            # this is an assumption about what is going on here
            return "insertion", new_coord
    return "none", 0


def error_adjusting(target, consensus, donor_VSG, donor_start, donor_end, error, donor_VSGs):
    """
    This function runs multiple versions of the seq_comparison to determine which works the best
    and returns the fewest identified recombination events. - representing the best alignment
    Args:
        target: a string, the target sequence
        consensus: a string, the consensus sequence
        donor_VSG: a string, the VSG donor name
        donor_start: an int, the donor VSG start position
        donor_end: an int, the donor VSG end position
        error: a string, a description of the error
        donor_VSGs: a dictionary, all the donor VSGs. key = VSG name, value = sequence

    Returns:
        cross_overs: a list, all the cross over events identified
    """
    # based on the error, the direction of adjustment
    if error == "insertion":
        adjust = 1
    elif error == "deletion":
        adjust = -1

    # test all possible combinations of shifted alignments, the one with the fewest recombinations
    # is returned
    cross_overs_reg = seq_comparison(target, consensus,
                                     donor_VSGs[donor_VSG][
                                     donor_start:donor_end])
    cross_overs_longfront = seq_comparison(target, consensus,
                                           donor_VSGs[donor_VSG][
                                           donor_start + adjust:donor_end])
    cross_overs_longend = seq_comparison(target, consensus,
                                         donor_VSGs[donor_VSG][
                                         donor_start:donor_end + adjust])

    options = [cross_overs_reg, cross_overs_longfront, cross_overs_longend]
    len_errors = list(map(len, options))
    min_co_errors = len_errors.index(min(len_errors))
    cross_overs = options[min_co_errors]

    return cross_overs


def overlap_quality_check(part_seq1, part_seq2, qual_seq1, qual_seq2):
    """
    takes both sequences from overlaps and
    Args:
        part_seq1: sequence 5'-3' in the overlap
        part_seq2: sequence 5'-3' in the overlap
        qual_seq1: corresponding seq1 quality
        qual_seq2: corresponding seq2 quality

    Returns:
        combined_seq: sequence 5'-3' with seq chosen by quality score
        error_flag: a flag, False if no error detected, True if error detected
    """
    error_flag = False
    final_seq = list(part_seq1)
    qual_offset = 33
    edits = lev.editops(part_seq1, part_seq2)
    for base in edits:
        if "insert" in base[0]:
            return "", True
        elif "delete" in base[0]:
            return "", True
        qual1 = ord(qual_seq1[base[1]]) - qual_offset
        qual2 = ord(qual_seq2[base[1]]) - qual_offset
        if qual1 > qual2:
            final_seq[base[1]] = part_seq1[base[1]]
        else:
            final_seq[base[1]] = part_seq2[base[1]]
    final_seq_str = "".join(final_seq)
    return final_seq_str, error_flag


def read_in_consensus(file_name):
    """
    A function for reading in the consensus sequences created in define_read_consensus
    Args:
        file_name: the name of file in .fa format containing all the consensus sequences.

    Returns:
        sequences: a list, the tuples of the base with the read info and the consensus sequence

    """
    sequences = []
    with open(file_name, "r") as con_file:
        for line in con_file:
            if ">" in line:
                title = line
                sample_info = title.strip().split("_")
                mouse = sample_info[1]
                day = sample_info[2]
                genotype = sample_info[3]
                primer = sample_info[4]
                UMI = sample_info[5]
                count = sample_info[6].split(": ")[0]
                vsgs = sample_info[6].split(": ")[1].split("; ")[0]
                un_anch_seq = sample_info[6].split(": ")[1].split("; ")[1]
                read2_cutoff = sample_info[6].split(": ")[1].split("; ")[2]
                base = mouse + "_" + day + "_" + genotype + "_" + primer + "_" + UMI + "_" + count + \
                    "_" + vsgs + "_" + un_anch_seq + "_" + read2_cutoff
            else:
                con_seq = line
                sequences.append((base, con_seq))

    return sequences


def mosaic_find(primer_file, donor_VSG_file):
    """
    run through the consensus sequences and solo R1s to identify the recombination sites between
    the target and the donors. If multiple donors were possible, all donors are tested. If one turns
    out to be a best fit, only the best fit is reported. Otherwise multiple are reported.
    Args:
        primer_file: a string, the text file with the primers ending in .txt
            antat_primers.txt
        donor_VSG_file: a string, path to the FASTA file containing the donor VSGs

    Returns:

    """
    # variable initialization
    donor_VSGs = read_in_donorVSGs(donor_VSG_file)
    primer_dict, primer_seq_dict = aux_functions.read_in_primers(primer_file)
    target_primer_range = target_primer_id(primer_dict)
    output = "double_crosses.csv"
    output2 = "ind_cross_over_events.csv"
    consensus_seqs = "consensus.fa"
    r1_only_seqs = "no_consensus_R1.fa"
    summary_results = "output.txt"
    weird_error = 0
    target_counts = 0
    total = 0

    with open(output, "w") as double_crosses, open(output2, "w") as cross_over_output:

        # read in consensus reads and raw R1s
        consensus_list = read_in_consensus(consensus_seqs)
        r1_files = read_in_consensus(r1_only_seqs)
        consensus_list = consensus_list + r1_files

        # loop through each read and identify cross overs
        for read_entry in consensus_list:
            doubles_cross_overs = []
            cross_overs = []
            total += 1
            print(read_entry)

            # identify the read details
            if ", " in read_entry[0]:
                interim = read_entry[0].split("(")[1].split(")")[0].replace(", ", "")
                info = (read_entry[0].split("(")[0] + interim + read_entry[0].split(")")[1]).strip().split("_")
            else:
                info = read_entry[0].strip().split("_")
            mouse = info[0]
            day = info[1]
            genotype = info[2]
            primer = info[3]
            direction = info[3][1]
            UMI = info[4]
            count = info[5]
            un_anch_seq = info[7]
            con_seq = read_entry[1].strip()

            # isolate the candidates
            candidate_VSGs = info[6].replace("(", "").replace(")", "").split("Tbb1125VSG-")
            candidate_VSGs.pop(0)

            # test consensus against the target just in case errors taken care of through consensus
            allowed_mismatch = str(1)
            target_regex = "(?b)(" + con_seq + "){e<=" + allowed_mismatch + "}"
            target_match = regex.search(target_regex, gbl.full_length_target)

            if target_match:
                print("This is actually the target")
                target_counts += 1
                continue

            # allows us to dynamically remove donors found to not be a match
            cand_copy = candidate_VSGs.copy()

            multiple_donors = []
            for vsg in candidate_VSGs:
                # rename the VSG
                VSG_name = "VSG-" + vsg

                # locate the positions within the donor where the read can be found
                donor_start, donor_end, indel_error = get_donor_range(con_seq, un_anch_seq, VSG_name, donor_VSGs)

                # if the donor is a bad match after consensus building,
                # 0 and 0 are returned as donor start/end
                if donor_start == 0 and donor_end == 0:
                    cand_copy.remove(vsg)
                    continue

                # check for an error near the primer
                # based on errors, redefine the donor start location and end location
                early_error, new_consensus_pos = primer_gap_test(con_seq, primer)
                if early_error == "insertion":
                    print("insertion_detected")
                    if "F" in primer:
                        primer_length_est = 20
                        con_seq = con_seq[primer_length_est:]
                        donor_start = donor_start + primer_length_est
                    else:
                        primer_length_est = -20
                        con_seq = con_seq[:primer_length_est]
                        donor_end = donor_end + primer_length_est
                elif early_error == "deletion":
                    print("deletion_detected")
                    if "F" in primer:
                        primer_length_est = 20
                        con_seq = con_seq[primer_length_est:]
                        donor_start = donor_start + primer_length_est
                    else:
                        primer_length_est = -20
                        con_seq = con_seq[:primer_length_est]
                        donor_end = donor_end + primer_length_est

                # locate the positions within the target where the read can be found
                if "F" in primer:
                    target_start = target_primer_range[primer][0]
                    target_end = target_primer_range[primer][0] + len(con_seq)
                else:
                    target_start = target_primer_range[primer][1] - len(con_seq)
                    target_end = target_primer_range[primer][1]

                # based on the ranges determined, obtain the target and donor sequences
                target_pred_seq = gbl.full_length_target[target_start: target_end]
                donor_seq = donor_VSGs[VSG_name][donor_start:donor_end]

                # determine recombination sites if errors have been detected or if errors are absent
                if indel_error != "none":
                    cross_overs = error_adjusting(target_pred_seq, con_seq, VSG_name, donor_start, donor_end,
                                                  indel_error, donor_VSGs)
                else:
                    cross_overs = seq_comparison(target_pred_seq, con_seq, donor_seq)

                if len(con_seq) == len(target_pred_seq) == len(donor_seq):
                    pass
                else:
                    print(len(con_seq))
                    print(len(target_pred_seq))
                    print(len(donor_seq))
                    print(count)

                if len(cross_overs) > 4:
                    print("some error is here, read removed")
                    weird_error += 1

                if len(candidate_VSGs) >= 2:
                    multiple_donors.append((vsg, cross_overs, target_start))
            candidate_VSGs = cand_copy.copy()
            # if no donor identified with good recombinations associated with it, remove read
            if len(candidate_VSGs) == 0:
                continue
            # check for a preferred donor match, if multiple donors originally identified
            if len(multiple_donors) != 0:
                updated_list = narrow_down_donor(multiple_donors)
            else:
                updated_list = [(vsg, cross_overs, target_start)]

            # these are independent of what the donor is
            for entry in updated_list:
                doubles_cross_overs.append(identify_double_COs(cross_overs, target_pred_seq, con_seq,
                                                               target_start, entry[0]))
            if len(doubles_cross_overs[0]) != 0:
                # previously had len of individual_cross_overs not equal to zero, incase this is impt.
                for donor in doubles_cross_overs:
                    for double_cross_over in donor:
                        donor_vsg = "Tbb1125VSG-" + str(double_cross_over[7])
                        double_crosses.write(", ".join([mouse, day, genotype,  primer, direction,
                                                            count, str(double_cross_over[0]),
                                                            str(double_cross_over[1][2]), str(double_cross_over[1][0]),
                                                            str(double_cross_over[1][1]), str(double_cross_over[2][2]),
                                                            str(double_cross_over[2][0]), str(double_cross_over[2][1]),
                                                            str(double_cross_over[3]), str(double_cross_over[4]),
                                                            str(double_cross_over[5]), str(double_cross_over[6]),
                                                            str(donor_vsg)]) + "\n")
            if updated_list[0][2] >= 0:
                for mosaic in updated_list:
                    vsg_donor = "Tbb1125VSG-" + str(mosaic[0])
                    for co_event in mosaic[1]:
                        cross_over_output.write(", ".join([mouse, day, genotype, primer, direction,
                                                               count, vsg_donor, str(mosaic[2]),
                                                               str(co_event[0]), str(co_event[1]),
                                                               co_event[2], UMI]) + "\n")
            candidate_VSGs = cand_copy.copy()
            if len(candidate_VSGs) == 0:
                continue


def main():
    """Execute the functions"""

    mosaic_find("antat_primers.txt", "EATRO1125_vsgs_long.fa")


if __name__ == '__main__':
    main()

