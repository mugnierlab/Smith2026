"""AMP Seq: This module creates mosaic consensus sequences from overlapping reads.
   Smith 2024
   Author: Jaclyn Smith """
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import aux_functions
import Levenshtein as lev
import regex
import glob
import global_target as gbl
import vsg_align_supplement


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
            cor_seq = vsg_align_supplement.rev_complement(primer_seq)
            primer_regex = "(?b)(" + cor_seq + "){e<=" + allowed_mismatch + "}"
            primer_return = regex.search(primer_regex, gbl.full_length_target)
            end = primer_return.span()[1]
            start = end - 150

        primer_ranges[primer_name] = (start, end)
    return primer_ranges


def get_consensus(seq_1, seq_2, qual_1, qual_2, un_anch_seq, primer, target_primer_range):
    """ create consensus sequence from both reads
        This function determines if there is overlap between the two reads and then combines the reads into
        one consensus sequence
        Args:
            seq_1: string, sequence of read1
            seq_2: string, sequence of read2
            qual_1: string, the quality  scores of read1, in 5'-3'
            qual_2: string, the quality scores of read2, in 5'-3'
            un_anch_seq: string, sequence of the donorVSG part of read1, sense with target
            primer: string, the primer used to amplify the fragment
            target_primer_range: a dictionary, key = primer name; value = tuple with primer range

        Returns:
            read_consensus: string, consensus sequence of combined R1 and R2, returns "" if no consensus exists
            consensus_donor_seq: string, consensus donor VSG sequence, returns un_anch_seq if no consensus exists
        """

    frags = seq_1.split(un_anch_seq)
    allowed_mismatch = str(1)
    recheck = False

    if frags[1] == "" and frags[0] == "":
        # None of the target was trimmed from either end, as a result, a consensus sequence is difficult to generate
        # or may be impossible due to spacing differences between the two sequenced fragments
        read_consensus = ""
        recheck = True
        return read_consensus, recheck

    elif frags[1] == "" and frags[0] != "":
        read1_target = frags[0]
    elif frags[0] == "" and frags[1] != "":
        read1_target = frags[1]
    elif frags[0] != "" and frags[1] != "":
        if len(frags[1]) >= len(frags[0]):
            read1_target = frags[1]
        else:
            read1_target = frags[0]

    if len(read1_target) < 20:
        # I think this will not yield anything by default, but just in case -
        # I do not want shorter than 20 just to make sure that it gets the right spot.
        print("too short for consensus")
        read_consensus = ""
        recheck = True
        return read_consensus, recheck

    target_return_regex = "(?b)(" + read1_target + "){e<=" + allowed_mismatch + "}"
    target_return = regex.search(target_return_regex, gbl.full_length_target)
    if target_return:
        target_return_range = target_return.span()
    else:
        print("error in trim, bad read")
        read_consensus = ""
        recheck = True
        return read_consensus, recheck

    # check if reads overlap post location identification
    r2 = range(target_primer_range[primer][0], target_primer_range[primer][1])
    r1 = range(target_return_range[0], target_return_range[1])
    read_set = set(r2)
    if len(read_set.intersection(r1)) == 0:
        print("these reads don't overlap")
        read_consensus = ""
        return read_consensus, recheck

    if read1_target == frags[1]:
        read_1_end = target_return_range[1]
    else:
        read_1_end = target_return_range[0] + len(seq_1)

    if "F" in primer:
        print(seq_1)
        read_consensus, error = consensus_fwd(seq_1, seq_2, qual_1, qual_2, read_1_end, primer, target_primer_range)
    else:
        read_consensus, error = consensus_rev(seq_1, seq_2, qual_1, qual_2, read_1_end, primer, target_primer_range)

    if error:
        recheck = True
    return read_consensus, recheck


def consensus_fwd(seq1, seq2, qual1, qual2, read1_end, primer, target_primer_range):
    """ create consensus sequence from both reads
    This function determines if there is overlap between the two reads and then combines the reads into
    one consensus sequence
    Args:
        seq1: string, sequence of read1
        seq2: string, sequence of read2
        qual1: string, the quality  scores of read1, in 5'-3'
        qual2: string, the quality scores of read2, in 5'-3'
        read1_end: int, position, the end of the read1 area
        primer: string, primer used to amplify fragment
        target_primer_range: a dictionary, key = primer name; value = tuple with primer range

    Returns:
        read_consensus: string, consensus sequence of combined R1 and R2, returns "" if no consensus found
    """
    error_detected = False
    read_consensus = seq2[:read1_end - target_primer_range[primer][0] - len(seq1)]
    read1_trim = len(seq2) - (read1_end - target_primer_range[primer][0] - len(seq1))
    # then add on overlapping piece
    part1 = seq1[:read1_trim]
    part2 = seq2.upper()[read1_end - target_primer_range[primer][0] - len(seq1):]
    part1q = qual1[:read1_trim]
    part2q = qual2[read1_end - target_primer_range[primer][0] - len(seq1):]
    if len(part1) != len(part2):
        # this is probably from some kind of improper trimming
        part2 = part2[:len(part1)]
    con_overlap, error_flag = overlap_quality_check(part1, part2, part1q, part2q)
    if error_flag:
        print("some unknown error occurred. Previous has been a deletion")
        read_consensus = ""
        error_detected = True
        return read_consensus, error_detected
    read_consensus += con_overlap

    # finally add on seq1 only piece
    read_consensus += seq1[read1_trim:]

    return read_consensus, error_detected


def consensus_rev(seq1, seq2, qual1, qual2, read1_end, primer, target_primer_range):
    """ create consensus sequence from both reads
    This function determines if there is overlap between the two reads and then combines the reads into
    one consensus sequence
    Args:
        seq1: string, sequence of read1
        seq2: string, sequence of read2
        qual1: string, the quality  scores of read1, in 5'-3'
        qual2: string, the quality scores of read2, in 5'-3'
        read1_end: int, position, the end of the read1 area
        primer: string, the primer used to amplify the fragment
        target_primer_range: a dictionary, key = primer name; value = tuple with primer range

    Returns:
        read_consensus: string, consensus sequence of combined R1 and R2, returns "" if no consensus found
    """
    error_detected = False
    allowed_mismatch = str(1)
    test_seq = seq2[len(seq2) - 20:]
    new_target_regex = "(?b)(" + test_seq + "){e<=" + allowed_mismatch + "}"
    new_target = regex.search(new_target_regex, gbl.full_length_target)
    if new_target:
        target_end_range = new_target.span()[1]
    else:
        read_consensus = ""
        error_detected = True
        return read_consensus, error_detected

    read_consensus = seq2[len(seq2) - (target_end_range - read1_end):]
    read1_trim = len(seq1) - (len(seq2) - len(read_consensus))

    part1 = seq1[read1_trim:]
    part2 = seq2.upper()[:len(seq2) - len(read_consensus)]
    part1q = qual1[read1_trim:]
    part2q = qual2[:len(seq2) - len(read_consensus)]

    if len(part1) != len(part2):
        part2 = part2[len(part2) - len(part1):]

    con_overlap, error_flag = overlap_quality_check(part1, part2, part1q, part2q)
    read_consensus = con_overlap + read_consensus

    read_consensus = seq1[:read1_trim] + read_consensus

    if error_flag:
        print("some unknown error occurred. Previously this has been a deletion!")
        read_consensus = ""
        error_detected = True
        return read_consensus, error_detected
    return read_consensus, error_detected


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


def unknown_overlap(seq1, seq2, qual1, qual2, primer):
    """
    This function tries to find overlaps between reads which are unrelated to target
     and create a consensus sequence
    Args:
        seq1: string, sequence of read1
        seq2: string, sequence of read2
        qual1: string, the quality  scores of read1, in 5'-3'
        qual2: string, the quality scores of read2, in 5'-3'
        primer: string, primer used to amplify fragment

    Returns:
        read_consensus: string, consensus sequence of combined R1 and R2, returns "" if no consensus found

    """
    query_len = 20
    allowed_mismatch = str(1)
    if "F" in primer:
        test_seq = seq2[len(seq2) - query_len:]
    else:
        test_seq = seq2[: query_len]

    target_return_regex = "(?b)(" + test_seq + "){e<=" + allowed_mismatch + "}"
    target_return = regex.search(target_return_regex, seq1)
    if target_return:
        target_return_range = target_return.span()
    else:
        print("these reads do not overlap.")
        read_consensus = ""
        return read_consensus

    if "F" in primer:
        read_1_split = target_return_range[1]
        # front part
        read_consensus = seq2[:len(seq2) - read_1_split]
        # overlap part
        part1 = seq1[:read_1_split]
        part2 = seq2[len(seq2) - read_1_split:]
        part1q = qual1[:read_1_split]
        part2q = qual2[len(seq2) - read_1_split:]

        overlap_seq, error = overlap_quality_check(part1, part2, part1q, part2q)

        if error:
            read_consensus = ""
            return read_consensus
        read_consensus += overlap_seq

        # back part
        read_consensus += seq1[read_1_split:]

    else:
        read_1_split = target_return_range[0]
        # front part
        read_consensus = seq1[:read_1_split]

        # overlap part
        part1 = seq1[read_1_split:]
        part2 = seq2[:len(seq1) - read_1_split]
        part1q = qual1[read_1_split:]
        part2q = qual2[:len(seq1) - read_1_split]

        overlap_seq, error = overlap_quality_check(part1, part2, part1q, part2q)

        if error:
            read_consensus = ""
            return read_consensus
        read_consensus += overlap_seq

        # back part
        read_consensus += seq2[len(seq1) - read_1_split:]
    return read_consensus


def consensus_builder():
    # get target primer ranges
    primer_dict, primer_seq_dict = aux_functions.read_in_primers("antat_primers.txt")
    target_ranges = target_primer_id(primer_dict)

    # output variable initialization
    count = 0
    mosaic_folder = "./mosaics/"
    con_R1 = "consensus_R1.fa"
    con_R2 = "consensus_R2.fa"
    output_consensus = "consensus.fa"
    no_con_R1 = "no_consensus_R1.fa"
    no_con_R2 = "no_consensus_R2.fa"
    files_to_check = "consensus_by_hand.fa"


    with open(con_R1, "w") as con_read1, open(con_R2, "w") as con_read2, \
            open(output_consensus, "w") as con_read, open(no_con_R1, "w") as no_read1, \
            open(no_con_R2, "w") as no_read2, open(files_to_check, "w") as by_hand:
        for consol_read in glob.glob(mosaic_folder + "*_reads_R1.fq"):
            print(consol_read)
            base = consol_read.split("/")[2].split("_mosaic_reads")[0]
            read1 = mosaic_folder + base + "_mosaic_reads_R1.fq"
            read2 = mosaic_folder + base + "_mosaic_reads_R2.fq"
            mosaic_possibilities = mosaic_folder + base + "_mosaic_info.txt"
            sample_info = base.split("_")
            primer = sample_info[3]
            direction = sample_info[3][1]

            with open(read2, "r") as seq_file_read2, open(read1, "r") as seq_file_read1, \
                    open(mosaic_possibilities, "r") as file3:
                for (title_1, sequence_1, quality_1), (title_2, sequence_2, quality_2), line in zip(
                        FastqGeneralIterator(seq_file_read1),
                        FastqGeneralIterator(seq_file_read2),
                        file3):
                    info = line.strip().split("; ")
                    UMI = title_1.strip().split(":")[-1][0:-2]
                    if direction == "F":
                        seq_1 = rev_complement(sequence_1).upper()
                        seq_2 = sequence_2.upper()
                        qual_1 = quality_1[::-1]
                        qual_2 = quality_2
                    else:
                        seq_1 = sequence_1.upper()
                        seq_2 = rev_complement(sequence_2).upper()
                        qual_1 = quality_1
                        qual_2 = quality_2[::-1]

                    candidate_VSGs = info[0].replace("(", "").replace(")", "").split("Tbb1125VSG-")
                    candidate_VSGs.pop(0)
                    consensus, recheck_alert = get_consensus(seq_1, seq_2, qual_1, qual_2, info[1], primer, target_ranges)
                    count += 1

                    if recheck_alert:
                        consensus = unknown_overlap(seq_1, seq_2, qual_1, qual_2, primer)
                        if consensus == "":
                            by_hand.write(">seq1_" + base + "_" + UMI + "_" + str(count) + ": " + line)
                            by_hand.write(sequence_1 + "\n")
                            by_hand.write(">seq2_" + base + "_" + UMI + "_" + str(count) + ": " + line)
                            by_hand.write(sequence_2 + "\n")

                    if consensus != "":
                        con_read1.write(">seq1_" + base + "_" + UMI + "_" + str(count) + ": " + line)
                        con_read1.write(sequence_1 + "\n")
                        con_read2.write(">seq2_" + base + "_" + UMI + "_" + str(count) + ": " + line)
                        con_read2.write(sequence_2 + "\n")
                        con_read.write(">consensus_" + base + "_" + UMI + "_" + str(count) + ": " + line)
                        con_read.write(consensus + "\n")
                    else:
                        no_read1.write(">seq1_" + base + "_" + UMI + "_" + str(count) + ": " + line)
                        no_read1.write(sequence_1 + "\n")
                        no_read2.write(">seq2_" + base + "_" + UMI + "_" + str(count) + ": " + line)
                        no_read2.write(sequence_2 + "\n")


def main():
    consensus_builder()


if __name__ == '__main__':
    main()
