"""VSG-AMP-Seq: This set of functions aligns reads to target VSG and identifies mosaic VSG reads candidates
   Smith 2024
   Author: Jaclyn Smith
"""
import glob
import os
import csv
from contextlib import ExitStack
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from progressbar import ProgressBar
import Levenshtein as lev
import regex
import vsg_align_supplement
import global_target as gbl


def mispriming(sequence_2, pred_seq):
    """
    Determines if this read is misassigned or primed in the wrong spot.
    Args:
        sequence_2: a string, sequence for read2
        pred_seq: a string, target predicted sequence which the primer should amplify

    Returns:
        r2_vsg_match, a regex object, if there is a match it is filled with match parameters
            if there is no match, None is returned
        retain, a boolean, returns True if read to be kept, but False if read not a match
        reason, a string, returns "short" or "misprimed"
    """
    retain = True
    seq_no_primer = sequence_2.translate(gbl.table)
    pred_seq_no_primer = pred_seq.translate(gbl.table)
    if len(seq_no_primer) <= 15:
        r2_vsg_match = None
        retain = False
        reason = "short"
        return r2_vsg_match, retain, reason
    if lev.distance(seq_no_primer[0:10], pred_seq_no_primer[0:10]) > 4:
        r2_vsg_match = None
        retain = False
        reason = "misprimed"
        return r2_vsg_match, retain, reason
    else:
        r2_n_num = seq_no_primer.count("N")
        allowed_mismatch = str(r2_n_num + 1)
        read_regex = "(?b)(" + seq_no_primer + "){s<=" + allowed_mismatch + "}"
        r2_vsg_match = regex.search(read_regex, gbl.full_length_target)
        reason = "none"
        return r2_vsg_match, retain, reason


def trimming_r1_check(sequence_2, un_anch_seq, primer):
    """
    A check for proper adapter trimming. If there trimming was not successful, a trimmed version of
    the read is returned. Checks first for overlap with R2. If not identified, also checks for the
    adapter sequence itself.
    Args:
        sequence_2: a string, the sequence of read2, but the 5' to 3' version of read2
        un_anch_seq: a string, the 5' to 3' version of read1, matching the direction of the target
        primer: a string, the direction of the primer used to amplify fragment
    Returns:
        un_anch_seq: a string, the updated version of read1, trimmed if necessary
        trimmed: boolean, True if trimming occurred, but False if kept intact.
    """
    # first, correct direction of sequence
    if "R" in primer:
        first_part = vsg_align_supplement.rev_complement(sequence_2.upper()[0:20])
    else:
        first_part = sequence_2.upper()[0:20]
    #first_part = seq_2.upper()[0:20]

    # search for first bps of seq2 in seq1
    read_regex = "(?b)(" + first_part + "){e<=" + "1" + "}"
    start_match = regex.search(read_regex, un_anch_seq)

    # if the first bps are within seq1, look further for proper adapter trimming
    if start_match:
        pos_trim = start_match.span()[0]
        un_anch_seq = un_anch_seq[pos_trim:len(un_anch_seq)]
        if pos_trim == 0:
            trimmed = False
        else:
            trimmed = True
    else:
        # if first bps not present, look instead for the adapter segments
        adapter_constant = "CCTCTCTATGGGCAGTCGGTGAT"
        read_regex = "(?b)(" + adapter_constant + "){e<=" + "2" + "}"
        constant_match = regex.search(read_regex, un_anch_seq)
        # if found, remove the adapter segments
        if constant_match:
            pos_trim = constant_match.span()[1]
            un_anch_seq = un_anch_seq[pos_trim:len(un_anch_seq)]
            trimmed = True
        else:
            trimmed = False

    return un_anch_seq, trimmed


def read_split(types_of_edits, sequence):
    """
    This function splits non-target reads into two parts. It takes into account Ns by allowing Ns
    to not count towards the error rate. The split occurs where up to one error is included in the
    target portion. If the read shifts from matching the target to having numerous discrepancies
    in a row (either errors or a new VSG), those are removed from the portion matching the target.
    This function only works in the direction of the tuple list, be careful with F/R fragments.
    Args:
        types_of_edits:list of tuples, this list, returned as fuzzy positions by regex gives a
            list of modifications to change one string to the other. in this case, our target
            to the read.
            pos 0 gives the type of change as a string.
            pos 1 gives position in target, pos 2 gives position in read.
        sequence: string, read sequence

    Returns:
        read_vsg_end: the position of the end of the target, start of the potential mosaic.
    """
    # initialize variables
    read_length = len(sequence)
    read_first_mut = True
    read_vsg_end = 0
    read_first_mut_pos = 0

    # run through edits. One mutation is allowed and still will count as target if not directly
    # followed by a second.
    for edit in types_of_edits:
        if read_length == edit[2]:
            read_vsg_end = read_length
            break
        # disregard Ns from counting as errors
        if sequence[edit[2]] != "N":
            # identify first mutation
            if read_first_mut:
                read_first_mut = False
                read_first_mut_pos = edit[1]
            else:
                # check if first mutation directly followed by a second.
                if read_first_mut_pos + 1 == edit[1]:
                    read_vsg_end = read_first_mut_pos
                else:
                    read_vsg_end = edit[2]
                break

        read_vsg_end = read_length

    return read_vsg_end


def read_split_reverse(types_of_edits, sequence):
    """
    This function splits non-target reads into two parts. It takes into account Ns by allowing
    Ns to not count towards the error rate. The split occurs where up to one error is included
    in the target portion. If the read shifts from matching the target to having numerous
    discrepancies in a row (either errors or a new VSG), those are removed from the portion
    matching the target. This function is for the opposite end of the strand as read_split.
    However, this function only works in the direction of the tuple list, be careful with
    F/R fragments.

    Args:
        types_of_edits:list of tuples, this list, returned as fuzzy positions by regex gives
            a list of modifications to change one string to the other. in this case, our target
            to the read.
            pos 0 gives the type of change as a string.
            pos 1 gives position in target, pos 2 gives position in read.
        sequence: string, read sequence

    Returns:
        read_vsg_end: the position of the end of the target, start of the potential mosaic.
    """
    # initialize variables
    read_length = len(sequence)
    read_first_mut = True
    read_vsg_end = 0
    read_first_mut_pos = 0

    # run through edits. One mutation is allowed and still will count as target if not directly
    # followed by a second.
    for edit in reversed(types_of_edits):
        if read_length == edit[2]:
            read_vsg_end = read_length
            break
        # disregard Ns from counting as errors
        if sequence[edit[2]] != "N":
            # identify first mutation
            if read_first_mut:
                read_first_mut = False
                read_first_mut_pos = edit[2]
            else:
                # check if first mutation directly followed by a second.
                if read_first_mut_pos - 1 == edit[2]:
                    read_vsg_end = read_first_mut_pos + 1

                else:
                    read_vsg_end = edit[2] + 1

                break
        read_vsg_end = read_length

    return read_vsg_end


def alternative_amplified_VSG(read1, primer, allowed_mismatch, vsgs_amped_by_primer):
    """
    This function searches through list of VSGs where the primer is likely to bind even if
    they are not combined with the target into a mosaic. They would be false positives
    that need to be removed. The first matching VSG in this case is returned.
    Args:
        read1: a string, un_anch_seq; R1 in 5' - 3' direction
        primer: a string, the primer used to amplify the fragment
        allowed_mismatch: an int, number of allowed mismatches to count as aligned
        vsgs_amped_by_primer: nested dictionary, key = primer, value = dictionary with
            key = VSG_name and value = sequence

    Returns:
        read_present: regex object, None if nothing is found, but a regex object if read is
            found within list
    """
    # look for reads aligning to alternative VSGs with primer sequence
    for title in vsgs_amped_by_primer[primer]:
        read_regex = "(?b)(" + read1 + "){e<=" + str(allowed_mismatch) + "}"
        read_present = regex.search(read_regex, vsgs_amped_by_primer[primer][title])
        if read_present:
            return read_present
    read_present = None
    return read_present


def vsg_mispriming_test(vsg, sequence_2, direction):
    """
    This function tests for vsgs that align to a putative donor bc they were probably misamplified
    off of that VSG rather than representing a mosaic between the donor and the target. There are
    cases where just the primer doesn't match up well, but all remaining sequences appear to be
    from another VSG. These should be removed.

    This is used after a positive identification of a putative donorVSG, to check R2 also doesn't
    completely match that VSG too.
    Args:
        vsg: a string, the VSG sequence
        sequence_2: a string, the read R2
        consol_read: a string, the title

    Returns:
        retain: a boolean, False if read is a mispriming event.
    """
    retain = True
    seq_no_primer = sequence_2.translate(gbl.table)
    if direction == "F":
        seq_2 = seq_no_primer
    else:
        seq_2 = vsg_align_supplement.rev_complement(seq_no_primer)

    allowed_mismatch = str(1)
    read_regex = "(?b)(" + seq_2 + "){e<=" + allowed_mismatch + "}"
    read_present = regex.search(read_regex, vsg)

    if read_present:
        retain = False
    return retain


def vsg_search(read, primer, allowed_mismatch, sequence_2, direction, vsg_seqs_by_primer):
    """
    Uses liberally trimmed un_anch_seq to find donorVSGs which match the inserted segment
    into the target. Returns a list of possible candidates.
    Args:
        read: a string, un_anch_seq; R1 in 5'-3'
        primer: a string, primer used to obtain sequence
        allowed_mismatch: an int, number of mismatches permitted before a sequence does not align.
        sequence_2: a string, R2
        direction: a string, the direction of primer used to amplify the fragment
        vsg_seqs_by_primer: nested dictionary, key = primer, value = dictionary with
            key = VSG_name and value = sequence

    Returns:
        read_in_vsgs: a list, a list of candidate donor VSGs
        likely_misprimed: a boolean, True if misprimed
    """
    # look for reads aligning to alternative VSGs with primer sequence
    read_in_vsgs = []
    likely_misprimed = False
    for title in vsg_seqs_by_primer[primer]:
        read_regex = "(?b)(" + read + "){e<=" + str(allowed_mismatch) + "}"
        read_present = regex.search(read_regex, vsg_seqs_by_primer[primer][title])
        if read_present:
            keep_read = vsg_mispriming_test(vsg_seqs_by_primer[primer][title],
                                            sequence_2, direction)
            if keep_read:
                read_in_vsgs.append(title)
            else:
                likely_misprimed = True
    return read_in_vsgs, likely_misprimed


def anchored_end_trim(sequence_2, un_anch_seq, r2_vsg_end, r2_distance, direction):
    """
    This method trims un_anch_seq (5'-3' R1) at the anchored end. A trim position is determined
    based on 1) partial overlap with R2 known to be the target and 2) position based on a sequence
    within un_anch_seq. Each of these is done independently, then the  more liberal trim is chosen.
    Once a known position is determined, direction comparison determined which part of the read
    matches the target VSG and can be removed. This method works for both F/R reads.
    Trimmed product is checked for also matching to target, if so, no trim is performed and
    double_vsg will return as True.

    Args:
        sequence_2: a string, DNA sequence of sequence 2 in either strand depending on primer
            direction
        un_anch_seq: a string, the part of sequence 1 in 5' to 3' direction which will be trimmed.
        r2_vsg_end: an int, a position in seq_2 where target ends
        r2_distance: an int, number of differences between read2 and target
        direction: a string: either "R" or "F", the direction of the primer used to obtain the
            fragment

    Returns:
        vsg_trimmed_seq: a string, the DNA sequence of the final trimmed version of un_anch_seq
        double_vsg: a boolean: True if both fragments after trim align to target VSG
    """
    # initialize variables
    cut = False
    double_vsg = False
    read_overlap_len = 10
    overlap = False
    search_len = 20
    allowed_mismatch = str(1)
    if direction == "F":
        seq2_r1_vsg_end = 0
        loc_r1_vsg_end = 0
    else:
        seq2_r1_vsg_end = 150
        loc_r1_vsg_end = 150

    # predict trimming via overlap
    last_part = sequence_2[r2_vsg_end - read_overlap_len: r2_vsg_end]
    if direction == "R":
        last_part = vsg_align_supplement.rev_complement(last_part)

    read_regex = "(?b)(" + last_part + "){e<=" + allowed_mismatch + "}"
    overlap_match = regex.search(read_regex, un_anch_seq)
    if len(last_part) == 0:
        overlap_match = None

    if overlap_match:
        overlap = True
        if direction == "F":
            seq2_r1_vsg_end = overlap_match.span()[1]
        else:
            seq2_r1_vsg_end = overlap_match.span()[0]

    if direction == "F":
        search_seq = un_anch_seq[:search_len]
    else:
        search_seq = un_anch_seq[len(un_anch_seq) - search_len:]

    # predict  trimming via position
    read_regex = "(?b)(" + search_seq + "){e<=" + allowed_mismatch + "}"
    location_match = regex.search(read_regex, gbl.full_length_target)
    if location_match:
        if direction == "F":
            loc_r1_vsg = location_match.span()[0]
            r1_frag_end = loc_r1_vsg + len(un_anch_seq)
        else:
            r1_frag_end = location_match.span()[1]
            loc_r1_vsg = r1_frag_end - len(un_anch_seq)

        pred_vsg_seq = gbl.full_length_target[loc_r1_vsg:r1_frag_end]
        pred_vsg_r1_edits = lev.editops(pred_vsg_seq, un_anch_seq)
        if direction == "F":
            loc_r1_vsg_end = read_split(pred_vsg_r1_edits, un_anch_seq)
        else:
            loc_r1_vsg_end = read_split_reverse(pred_vsg_r1_edits, un_anch_seq)

    # compare methods and choose more liberal trim
    # allowed_mismatch = str(2)
    if direction == "F":
        r1_vsg_end = max(seq2_r1_vsg_end, loc_r1_vsg_end)
        frag_vsg = un_anch_seq[:r1_vsg_end]
        frag_other = un_anch_seq[r1_vsg_end:]
    else:
        r1_vsg_end = min(seq2_r1_vsg_end, loc_r1_vsg_end)
        frag_vsg = un_anch_seq[r1_vsg_end:]
        frag_other = un_anch_seq[:r1_vsg_end]

    # confirm trim should be executed. Check for both vsg part matching target and
    # non-vsg frag not matching target
    # only make a cut if it looks good.
    vsg_frag_regex = "(?b)(" + frag_vsg + "){e<=" + allowed_mismatch + "}"
    other_frag_regex = "(?b)(" + frag_other + "){e<=" + allowed_mismatch + "}"
    vsg_frag_check = regex.search(vsg_frag_regex, gbl.full_length_target)
    # if R2 matches target VSG throughout the overlapping portion, extra observed errors
    # are likely from sequencing & can be ignored
    if vsg_frag_check or (overlap and r2_distance <= 1):
        other_frag_check = regex.search(other_frag_regex, gbl.full_length_target)
        if other_frag_check is None:
            cut = True
        else:
            double_vsg = True
    if cut:
        if direction == "F":
            vsg_trimmed_seq = un_anch_seq[r1_vsg_end:]
        else:
            vsg_trimmed_seq = un_anch_seq[:r1_vsg_end]
    else:
        vsg_trimmed_seq = un_anch_seq
    return vsg_trimmed_seq, double_vsg


def unanchored_end_trim(un_anch_seq, direction):
    """
    This method trims un_anch_seq (5'-3' R1) at the unanchored end. A trim position is determined
    by position based on a sequence within un_anch_seq. Once a known position is determined,
    direction comparison determined which part of the read matches the target VSG and can be
    removed. This method works for both F/R reads. Trimmed product is checked for also matching
    to target, if so, no trim is performed and double_vsg will return as True.
    Args:
        un_anch_seq: a string, the part of sequence 1 in 5' to 3' direction which will be trimmed.
    Returns:
        vsg_trimmed_seq: a string, the DNA sequence of the final trimmed version of un_anch_seq
        double_vsg: a boolean: True if both fragments after trim align to target VSG
    """
    cut = False
    double_vsg = False
    allowed_mismatch = str(1)
    test_length = 10
    vsg_test_len = 20

    if len(un_anch_seq) < test_length:
        pass
    else:
        if direction == "F":
            test_seq = un_anch_seq[len(un_anch_seq) - test_length:]
        else:
            test_seq = un_anch_seq[:test_length]
        read_regex = "(?b)(" + test_seq + "){e<=" + allowed_mismatch + "}"
        r1_vsg_far_end = regex.search(read_regex, gbl.full_length_target)
        if r1_vsg_far_end:
            if direction == "F":
                r1_vsg_end = r1_vsg_far_end.span()[1]
                r1_frag_start = r1_vsg_end - len(un_anch_seq)
            else:
                r1_frag_start = r1_vsg_far_end.span()[0]
                r1_vsg_end = r1_frag_start + len(un_anch_seq)
            pred_vsg_seq = gbl.full_length_target[r1_frag_start:r1_vsg_end]
            pred_vsg_r1_edits = lev.editops(pred_vsg_seq, un_anch_seq)

            if direction == "F":
                r1_trimmed_end = read_split_reverse(pred_vsg_r1_edits, un_anch_seq)
                frag_other = un_anch_seq[:r1_trimmed_end]
                frag_vsg_3 = un_anch_seq[r1_trimmed_end:]
            else:
                r1_trimmed_end = read_split(pred_vsg_r1_edits, un_anch_seq)
                frag_other = un_anch_seq[r1_trimmed_end:]
                frag_vsg_3 = un_anch_seq[:r1_trimmed_end]

            vsg_frag_regex = "(?b)(" + frag_vsg_3 + "){e<=" + allowed_mismatch + "}"
            if len(frag_vsg_3) >= vsg_test_len:  # 20
                vsg_frag_check = regex.search(vsg_frag_regex, gbl.full_length_target)
            else:
                vsg_frag_check = None
            other_frag_regex = "(?b)(" + frag_other + "){e<=" + allowed_mismatch + "}"
            if vsg_frag_check:
                if len(frag_other) >= vsg_test_len:  # 20
                    other_frag_check = regex.search(other_frag_regex, gbl.full_length_target)
                else:
                    other_frag_check = None
                if other_frag_check is None:
                    cut = True
                else:
                    double_vsg = True
    if cut:
        if direction == "F":
            vsg_trimmed_seq = un_anch_seq[:r1_trimmed_end]
        else:
            vsg_trimmed_seq = un_anch_seq[r1_trimmed_end:]
    else:
        vsg_trimmed_seq = un_anch_seq
    return vsg_trimmed_seq, double_vsg


def alt_splicing_remover(un_anch_seq, alt_splice_counts):
    """
    This function checks all identified versions of alt splicing as specified in global_target.py. Alt splicing may
    still be something of interest, but removing these reads from the pipeline is important for speed. Because of that,
    matching reads are written to a separate output file, but are counted as matching the target. This only applies to
    reverse reads.
    Args:
        un_anch_seq: a string, the part of sequence 1 in 5' to 3' direction
        alt_splice_counts: a dictionary, key = name ie splice_one; value = counts of sequences that correspond to the
            alt splicing sequence

    Returns:
        alt_found: a boolean, False if no alt splicing is found, True if alt splicing match is found
        alt_splice_counts: a dictionary, updated from arguments

    """
    allowed_mismatch = str(1)
    alt_found = False
    for entry in gbl.alt_splicing:
        read_regex = "(?b)(" + un_anch_seq + "){e<=" + allowed_mismatch + "}"
        r1_alt_match = regex.search(read_regex, gbl.alt_splicing[entry])

        # check for match
        if r1_alt_match:
            alt_splice_counts[entry] += 1
            alt_found = True
            break
    return alt_found, alt_splice_counts


def vsg_align(source_folder):
    """
    Aligns all reads to the target VSG to identify mosaic hybrid reads.
    All results are written to files.
    Args:
        source_folder: a string, the folder location for the reads to run this. 
    """

    vsg_seqs_by_primer = {}
    vsgs_amped_by_primer = {}

    for primer in gbl.primer_dict.values():
        vsg_seqs_by_primer[primer] = vsg_align_supplement.read_in_candidateVSGs(primer, "absent")
        vsgs_amped_by_primer[primer] = vsg_align_supplement.read_in_candidateVSGs(primer, "contains")

    alt_splice = {}
    for alt_vsg in gbl.alt_splicing:
        alt_splice[alt_vsg] = 0

    folder_target = "target_reads"
    try:
        os.mkdir(folder_target)
        print("folder '{}' created ".format(folder_target))
    except FileExistsError:
        print("folder {} already exists".format(folder_target))
    folder_target = "./" + folder_target + "/"

    folder_mosaics = "mosaics"
    try:
        os.mkdir(folder_mosaics)
        print("folder '{}' created ".format(folder_mosaics))
    except FileExistsError:
        print("folder {} already exists".format(folder_mosaics))
    folder_mosaics = "./" + folder_mosaics + "/"

    folder_unknowns = "mosaic_unknowns"
    try:
        os.mkdir(folder_unknowns)
        print("folder '{}' created ".format(folder_unknowns))
    except FileExistsError:
        print("folder {} already exists".format(folder_unknowns))
    folder_unknowns = "./" + folder_unknowns + "/"

    folder_invert = "inversion_events"
    try:
        os.mkdir(folder_invert)
        print("folder '{}' created ".format(folder_invert))
    except FileExistsError:
        print("folder {} already exists".format(folder_invert))
    folder_invert = "./" + folder_invert + "/"

    folder_double = "double_events"
    try:
        os.mkdir(folder_double)
        print("folder '{}' created ".format(folder_double))
    except FileExistsError:
        print("folder {} already exists".format(folder_double))
    folder_double = "./" + folder_double + "/"

    folder_alt_splice = "alt_splice"
    try:
        os.mkdir(folder_alt_splice)
        print("folder '{}' created ".format(folder_alt_splice))
    except FileExistsError:
        print("folder {} already exists".format(folder_alt_splice))
    folder_alt_splice = "./" + folder_alt_splice + "/"

    folder_alt_vsg = "alt_vsg"
    try:
        os.mkdir(folder_alt_vsg)
        print("folder '{}' created ".format(folder_alt_vsg))
    except FileExistsError:
        print("folder {} already exists".format(folder_alt_vsg))
    folder_alt_vsg = "./" + folder_alt_vsg + "/"

    folder_alt_prime = "alt_prime"
    try:
        os.mkdir(folder_alt_prime)
        print("folder '{}' created ".format(folder_alt_prime))
    except FileExistsError:
        print("folder {} already exists".format(folder_alt_prime))
    folder_alt_prime = "./" + folder_alt_prime + "/"

    # initialize output file and write the stats.
    # Must be done this way as numerous times lost progress and did make note of
    # stats in a file in between sample runs
    alignment_stats = "target_alignment_stats.csv"
    header = ["mouse", "day", "genotype", "primer", "direction", "total_reads",
              "adapter_trimmed_reads", "removed_reads", "removed_reads:too_short", "removed_reads:too_many_Ns",
              "removed_reads:inversion", "removed_reads:double_target", "removed_reads:target_2errors",
              "removed_reads:removed_at_R2", "removed_reads:amplified_alt_target", "removed_reads:amplified_alt_vsg",
              "removed_reads:alt_splice", "target_R1", "target_R2", "unknown_R1", "unknown_R2", "mosaics_identified"]
    with open(alignment_stats, "a") as align:
        writer = csv.writer(align)
        writer.writerow(header)

    alt_splice_stats = "target_alt_splicing.csv"
    header_alt = ["mouse", "day", "genotype", "primer", "direction", "total_reads"]
    with open(alt_splice_stats, "a") as splice_out:
        writer = csv.writer(splice_out)
        writer.writerow(header_alt)

    pbar = ProgressBar()

    #for consol_read in pbar(glob.glob("./consol_reads/C*_1.fq")):
    for consol_read in pbar(glob.glob("./" + source_folder + "/*_1.fq")):
    #for consol_read in pbar(glob.glob("./demultiplexed_reads/*_1.fq")):

        read1 = consol_read
        read2 = consol_read.split("_1.fq")[0] + "_2.fq"

        # record sample details
        info = read2.strip().split("/")[2].split("_")
        mouse = info[0]
        day = info[1]
        genotype = info[2]
        primer = info[4]
        direction = info[4][1]

        # initalize output files
        mouse_info_for_output = mouse + "_" + day + "_" + genotype + \
                                "_" + primer + "_"
        mosaic_R1 = folder_mosaics + mouse_info_for_output + "mosaic_reads_R1.fq"
        mosaic_R2 = folder_mosaics + mouse_info_for_output + "mosaic_reads_R2.fq"
        mosaic_info = folder_mosaics + mouse_info_for_output + "mosaic_info.txt"
        target_R1 = folder_target + mouse_info_for_output + "target_reads_R1.fq"
        target_R2 = folder_target + mouse_info_for_output + "target_reads_R2.fq"
        unknown_R1 = folder_unknowns + mouse_info_for_output + "unknown_R1.fq"
        unknown_R2 = folder_unknowns + mouse_info_for_output + "unknown_R2.fq"
        invert_R1 = folder_invert + mouse_info_for_output + "inversion_examples_R1.fq"
        invert_R2 = folder_invert + mouse_info_for_output + "inversion_examples_R2.fq"
        double_R1 = folder_double + mouse_info_for_output + "double_vsg_examples_R1.fq"
        double_R2 = folder_double + mouse_info_for_output + "double_vsg_examples_R2.fq"
        alt_splice_R1 = folder_alt_splice + mouse_info_for_output + "alt_splice_R1.fq"
        alt_splice_R2 = folder_alt_splice + mouse_info_for_output + "alt_splice_R2.fq"
        alt_vsg_R1 = folder_alt_vsg + mouse_info_for_output + "alt_vsg_R1.fq"
        alt_vsg_R2 = folder_alt_vsg + mouse_info_for_output + "alt_vsg_R2.fq"
        alt_prime_R1 = folder_alt_prime + mouse_info_for_output + "alt_prime_R1.fq"
        alt_prime_R2 = folder_alt_prime + mouse_info_for_output + "alt_prime_R2.fq"

        # remove for testing
        if os.path.exists(mosaic_R1):
            continue

        # initialize counting variables
        total_reads = 0
        not_vsg_r2 = 0
        vsg_r2 = 0
        not_vsg_r1 = 0
        vsg_r1 = 0
        trimmed_count = 0
        alt_splice_count = 0

        # counts from removed reads
        too_short = 0
        amp_alt_target = 0
        amp_alt_vsg = 0
        too_many_Ns = 0
        invert = 0
        double_vsg_count = 0  # restricted to two fragments with shifted alignment
        vsg_2_errors = 0  # reads that align if two errors were permitted
        # probably not useful as mosaics anyway
        removed_at_R2 = 0  # these may be included in other categories
        removed_reads_count = 0  # too short/Ns, misprimed, large error
        vsgs_possible_mosaics = []
        alt_splice_found = False

        # because of so many possible destinations, need to use ExitStack() to prevent too many statically nested
        # blocks error. Using dictionary of positions for readability of code
        output_files = [mosaic_R1, mosaic_R2, mosaic_info, target_R1, target_R2, invert_R1,
                        invert_R2, double_R1, double_R2, unknown_R1, unknown_R2,
                        alt_splice_R1, alt_splice_R2, alt_vsg_R1, alt_vsg_R2, alt_prime_R1, alt_prime_R2]
        nums = {
            mosaic_R1: 0, mosaic_R2: 1, mosaic_info: 2, target_R1: 3, target_R2: 4, invert_R1: 5,
            invert_R2: 6, double_R1: 7, double_R2: 8, unknown_R1: 9, unknown_R2: 10,
            alt_splice_R1: 11, alt_splice_R2: 12, alt_vsg_R1: 13, alt_vsg_R2: 14, alt_prime_R1: 15, alt_prime_R2: 16
        }

        with ExitStack() as stack:
            files = [stack.enter_context(open(fname, "w")) for fname in output_files]

            with open(read2, "r") as seq_file_read2, open(read1, "r") as seq_file_read1:
                for (title_1, sequence_1, quality_1), (title_2, sequence_2, quality_2) in zip(
                        FastqGeneralIterator(seq_file_read1),
                        FastqGeneralIterator(seq_file_read2)):

                    total_reads += 1
                    sequence_2 = vsg_align_supplement.lowercase_primer(primer, sequence_2)

                    # begin Read2 processing
                    r2_length = len(sequence_2)
                    r2_distance = lev.distance(gbl.primer_seq_dict[primer][0:r2_length], sequence_2)
                    r2_n_num = sequence_2.count("N")
                    r2_distance = r2_distance - r2_n_num

                    if r2_n_num > 5:
                        too_many_Ns += 1
                        removed_reads_count += 1
                        continue

                    # Start with read2, the anchored read. This has a very predictable
                    # sequence in the target. because of this, the read is directly
                    # compared to the predicted sequence. If there is not a match
                    # the read is checked for mispriming. If it passes that test, ie the
                    # read is truely not antat, the first position where the read doesn't
                    # match the predicted sequence is located and reported as r2_vsg_end

                    if r2_distance > 1:
                        r2_types_of_edits = lev.editops(gbl.primer_seq_dict[primer][0:r2_length],
                                                        sequence_2)
                        pred_seq = gbl.primer_seq_dict[primer][0:r2_length]

                        r2_vsg_match, retain, reason = mispriming(sequence_2, pred_seq)
                        if not retain:
                            if reason == "short":
                                too_short += 1
                            elif reason == "misprimed":
                                amp_alt_target += 1
                                files[nums[alt_prime_R1]].write("@%s\n%s\n+\n%s\n" % (title_1, sequence_1, quality_1))
                                files[nums[alt_prime_R2]].write("@%s\n%s\n+\n%s\n" % (title_2, sequence_2, quality_2))
                            removed_at_R2 += 1
                            removed_reads_count += 1
                            continue
                        if r2_vsg_match:
                            vsg_r2 += 1
                            r2_distance = sum(list(r2_vsg_match.fuzzy_counts))
                            r2_vsg_end = len(sequence_2)
                        else:
                            not_vsg_r2 += 1
                            r2_vsg_end = read_split(r2_types_of_edits, sequence_2)
                    else:
                        r2_vsg_end = r2_length
                        vsg_r2 += 1

                    # read1 processing

                    # initialize variables
                    # test read1 against full length of target
                    if "F_" in consol_read:
                        un_anch_seq = vsg_align_supplement.rev_complement(sequence_1)
                    else:
                        un_anch_seq = sequence_1
                    r1_n_num = sequence_1.count("N")
                    if r1_n_num > 5:
                        too_many_Ns += 1
                        removed_reads_count += 1
                        continue
                    allowed_mismatch = str(r1_n_num + 1)

                    # this leverages regex and searches for the best match within the full length target
                    # if match found, a regex object is returned.
                    # if a match doesn't exist withing the parameters given, a None object is returned
                    # s for mismatch
                    # e for all errors, i for insertions and d for deletions
                    read_regex = "(?b)(" + un_anch_seq + "){e<=" + allowed_mismatch + "}"
                    r1_vsg_match = regex.search(read_regex, gbl.full_length_target)

                    # If not a match, check adapter trimming was successful!
                    if r1_vsg_match is None:
                        un_anch_seq, trimmed = trimming_r1_check(sequence_2, un_anch_seq, primer)
                        if trimmed:
                            if len(un_anch_seq) < 20:
                                # some reads appear to bind to 47 VSGs if the length is 15
                                # which obvs is not ideal
                                too_short += 1
                                removed_reads_count += 1
                                continue
                            # redo regex search with newly trimmed fragment
                            read_regex = "(?b)(" + un_anch_seq + "){e<=" + allowed_mismatch + "}"
                            r1_vsg_match = regex.search(read_regex, gbl.full_length_target)
                            trimmed_count += 1

                    # get range of additional information from read1
                    # check number and location of mismatches are the same
                    if r1_vsg_match:
                        vsg_r1 += 1
                        files[nums[target_R1]].write("@%s\n%s\n+\n%s\n" % (title_1, sequence_1, quality_1))
                        files[nums[target_R2]].write("@%s\n%s\n+\n%s\n" % (title_2, sequence_2, quality_2))
                    else:
                        # it can and needs to go twice, but not to allow 2 errors now.

                        un_anch_seq, double_vsg = anchored_end_trim(sequence_2, un_anch_seq, r2_vsg_end,
                                                                                r2_distance, direction)
                        if double_vsg:
                            # check that error is actually 2, rather than a shifted double alignment
                            double_mismatch = str(2)
                            read_regex_double = "(?b)(" + un_anch_seq + "){e<=" + double_mismatch + "}"
                            double_match = regex.search(read_regex_double, gbl.full_length_target)
                            # only shifted double alignments are written to double file.
                            if not double_match:
                                files[nums[double_R1]].write("@%s\n%s\n+\n%s\n" % (title_1, sequence_1, quality_1))
                                files[nums[double_R2]].write("@%s\n%s\n+\n%s\n" % (title_2, sequence_2, quality_2))
                                double_vsg_count += 1
                            else:
                                vsg_2_errors += 1
                            removed_reads_count += 1
                            continue

                        un_anch_seq, double_vsg = unanchored_end_trim(un_anch_seq, direction)

                        if double_vsg:
                            # check that error is not actually 3, rather than a shifted double alignment
                            double_mismatch = str(2)
                            read_regex_double = "(?b)(" + un_anch_seq + "){e<=" + double_mismatch + "}"
                            double_match = regex.search(read_regex_double, gbl.full_length_target)
                            # only shifted double alignments are written to double file.
                            if not double_match:
                                files[nums[double_R1]].write("@%s\n%s\n+\n%s\n" % (title_1, sequence_1, quality_1))
                                files[nums[double_R2]].write("@%s\n%s\n+\n%s\n" % (title_2, sequence_2, quality_2))
                                double_vsg_count += 1
                            else:
                                vsg_2_errors += 1
                            removed_reads_count += 1
                            continue

                        if len(un_anch_seq) < 20:
                            # some reads appear to bind to 47 VSGs if the length is 15,
                            # which obvs is not ideal
                            too_short += 1
                            removed_reads_count += 1
                            continue

                        if "R_" in consol_read:
                            alt_found, alt_splice = alt_splicing_remover(un_anch_seq, alt_splice)
                            if alt_found:
                                files[nums[alt_splice_R1]].write("@%s\n%s\n+\n%s\n" % (title_1, sequence_1, quality_1))
                                files[nums[alt_splice_R2]].write("@%s\n%s\n+\n%s\n" % (title_2, sequence_2, quality_2))
                                alt_splice_count += 1
                                alt_splice_found = True
                                continue

                        other_vsg = alternative_amplified_VSG(un_anch_seq,
                                                                    primer,
                                                                    allowed_mismatch,
                                                                    vsgs_amped_by_primer)
                        if other_vsg:
                            amp_alt_vsg += 1
                            removed_reads_count += 1
                            files[nums[alt_vsg_R1]].write("@%s\n%s\n+\n%s\n" % (title_1, sequence_1, quality_1))
                            files[nums[alt_vsg_R2]].write("@%s\n%s\n+\n%s\n" % (title_2, sequence_2, quality_2))
                            continue

                        invert_check = vsg_align_supplement.rev_complement(un_anch_seq)
                        read_regex = "(?b)(" + invert_check + "){e<=" + allowed_mismatch + "}"
                        invert_match = regex.search(read_regex, gbl.full_length_target)
                        if invert_match:
                            files[nums[invert_R1]].write("@%s\n%s\n+\n%s\n" % (title_1, sequence_1, quality_1))
                            files[nums[invert_R2]].write("@%s\n%s\n+\n%s\n" % (title_2, sequence_2, quality_2))
                            invert += 1
                            removed_reads_count += 1
                            continue

                        allowed_mismatch = str(1)

                        candidates, vsg_misprime = vsg_search(un_anch_seq,
                                                              primer,
                                                              allowed_mismatch,
                                                              sequence_2, direction,
                                                              vsg_seqs_by_primer)

                        if vsg_misprime:
                            amp_alt_target += 1
                            files[nums[alt_prime_R1]].write("@%s\n%s\n+\n%s\n" % (title_1, sequence_1, quality_1))
                            files[nums[alt_prime_R2]].write("@%s\n%s\n+\n%s\n" % (title_2, sequence_2, quality_2))
                            removed_reads_count += 1
                            continue

                        if len(candidates) >= 1:
                            files[nums[mosaic_R1]].write("@%s\n%s\n+\n%s\n" % (title_1, sequence_1, quality_1))
                            files[nums[mosaic_R2]].write("@%s\n%s\n+\n%s\n" % (title_2, sequence_2, quality_2))
                            cand_string = ""
                            for entry in candidates:
                                cand_string += entry

                            files[nums[mosaic_info]].write(("(" + cand_string + ")" + "; " + un_anch_seq + "; " +
                                                 str(r2_vsg_end) + "\n"))
                        else:

                            files[nums[unknown_R1]].write("@%s\n%s\n+\n%s\n" % (title_1, sequence_1, quality_1))
                            files[nums[unknown_R2]].write("@%s\n%s\n+\n%s\n" % (title_2, sequence_2, quality_2))
                            not_vsg_r1 += 1
                        vsgs_possible_mosaics.append(candidates)

                    num_mosaics = 0
                    for entry in vsgs_possible_mosaics:
                        if len(entry) >= 1:
                            num_mosaics += 1
        with open(alignment_stats, "a") as align_stats:
            writer = csv.writer(align_stats)
            writer.writerow([mouse, day, genotype, primer, direction, total_reads, trimmed_count,
                             removed_reads_count, too_short, too_many_Ns, invert, double_vsg_count,
                             vsg_2_errors, removed_at_R2, amp_alt_target, amp_alt_vsg, alt_splice_count, vsg_r1, vsg_r2,
                             not_vsg_r1, not_vsg_r2, num_mosaics])
        if alt_splice_found:
            with open(alt_splice_stats, "a") as splice_output:
                writer = csv.writer(splice_output)
                writer.writerow([mouse, day, genotype, primer, direction, total_reads])
                splice_output = []
                for entry in alt_splice:
                    splice_output.append(entry)
                    splice_output.append(alt_splice[entry])
                writer.writerow(splice_output)


def main():
    """Execute the functions"""
    vsg_align("demultiplexed_reads")


if __name__ == '__main__':
    main()
