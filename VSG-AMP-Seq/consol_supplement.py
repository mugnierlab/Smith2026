"""VSG-AMP-Seq: This set of functions is for consolidating reads from a group of reads
   Smith 2024
   Author: Jaclyn Smith
"""
import itertools
from more_itertools import locate
import Levenshtein as lev

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


def combine_reads(read_list, qual_list, min_qual):
    """ combine reads into a consensus seq

    This takes list of reads and quals and uses them to construct
    a consensus sequence. From lists of reads that are 98% similar,
    a consensus read will be the length of the longest read, the most popular
    base at every position, where each base is only counted if it has a sufficient
    quality score. If no base is a majority, the base is called as N, with a qual
    score of zero. This can only deal with one read a time, not both together.

    Args:
        read_list: an ordered list with all the reads to be consolidated
        qual_list: an ordered list with all quals corresponding to the reads,
            quals are in ascii format
        min_qual: a value, given by the user in the args of pipeline
    """
    # initialize string to build
    combined_read = ""
    combined_quals = ""
    qual_offset = 33

    # loop through each letter one at a time: for a given position, bases are in the reads tuple
    # and qualities for each of those bases are in the quals tuple; it will loop to the longest read
    # and for shorter reads, an X will be used to hold position.
    # have to change qual scores of X to 40 from 0, to keep them in the running for the consensus seq
    # this is from "!" to "I"
    for reads, quals in zip(itertools.zip_longest(*read_list, fillvalue="X"),
                            itertools.zip_longest(*qual_list, fillvalue="I")):

        # convert reads to sets, only one of each base is present
        base = set(reads)

        # bases which are the same for all reads have length 1
        if len(base) == 1:
            # concatenate to read
            combined_read = combined_read + list(base)[0]

            # obtain the max read quality, since bases are the same
            # loop only through unique quality scores
            unique_quality = set(quals)
            max_score = 0
            for score in unique_quality:
                # convert ascii score into numeric score
                numeric_score = ord(score) - qual_offset
                max_score = max(max_score, numeric_score)
            # concatenate to quality, convert numeric score back to ascii
            combined_quals = combined_quals + chr(max_score + qual_offset)
        else:
            # for bases where all reads are not the same
            # create a dict to count each bases
            base_dict = dict.fromkeys(base, 0)

            # loop through all reads, count only ones with high quality above
            # user specified min_quality
            i = 0
            while i < len(quals):
                # convert ascii score into numeric score
                numeric_score = ord(quals[i]) - qual_offset
                if numeric_score > min_qual:
                    base_dict[reads[i]] = base_dict[reads[i]] + 1
                i += 1

            # if all counts for bases are equal, then make the base an N,
            # this will catch all bases that are too low quality to count as well
            if len(set(base_dict.values())) == 1:
                consensus_base = "N"
            else:
                # assign the base that has the highest count
                consensus_base = max(base_dict, key=base_dict.get)
            combined_read = combined_read + consensus_base

            # look only at qualities for the base selected
            # this will prevent one high quality being used for a base to which it
            # does not correspond
            read_local = list(locate(reads, lambda x: x == consensus_base))
            max_score = 0
            for local in read_local:
                # convert ascii score into numeric score
                numeric_score = ord(quals[local]) - qual_offset
                max_score = max(max_score, numeric_score)

            combined_quals = combined_quals + chr(max_score + qual_offset)

    # check for trailing Ns/Xs and remove. if more than two, report multiple frags likely
    n_count = combined_read.count("N")
    x_count = combined_read.count("X")

    combined_read = combined_read.split("X")[0]
    combined_read = combined_read.strip("N")
    combined_quals = combined_quals[:len(combined_read)]

    return combined_read, combined_quals


def consensus_checker(consensus, read_list, cutoff):
    """
    This will compare the final consensus sequence to all reads to see if there are
    rogue reads which do not match the consensus. I want to get an idea of how many collisions/
    bad reads are in my consolidated data
    Args:
        consensus: string, DNA consensus sequence from the reads
        read_list: list of strings, DNA seqs which were used to consolidate
        cutoff: int, cutoff of similar

    Returns:
        consen: int, a count of the mismatching reads

    """
    # initialize variables
    consen_count = 0

    # loop through reads and compare them to consensus seq, if at least one hits the threshold,
    # count is returned
    for item in read_list:
        distance = lev.distance(consensus[:len(item)], item)
        if len(consensus) == 0:
            consen_count += 1
        elif distance/len(consensus)*100 > cutoff:
            consen_count += 1
    return consen_count


def consensus_sequence(reads_list, min_qual, consen_cutoff):
    """ create consensus sequence for binned reads

    This function takes in a list of read tuples and returns a consensus read along
    with a list of reads which are represented by this read and a count of reads which
    are a poor match to consensus. A read is considered a bad match when it differs by
    more than 10% over a range of equal length. This particular function breaks the reads
    into parts which can be fed into the combine_reads() function. Then puts them back
    together and returns output.
    Args:
        reads_list: list of read tuples which are to be consolidated
        min_qual: the minimum quality score for a base to count towards the consensus sequence
        consen_cutoff: the min difference allowed for reads to be considered similar to consensus
    Returns:
        updated_consol_read: read tuple representing the list (sequence, quality, name)
        name_record: list of all names collapsed into consol read
        bad_seq1_count: count of read1s which are more than 10% different from consensus
        bad_seq2_count: count of read2s which are more than 10% different from consensus
    """
    #initialize variables
    seq1 = []
    seq2 = []
    qual1 = []
    qual2 = []
    title = []

    # create lists for read1/2 and corresponding quality
    for read in reads_list:
        both_seqs = read[1].split("_")
        seq1.append(both_seqs[0])
        seq2.append(both_seqs[1])
        both_quals = read[2].split("_")
        qual1.append(both_quals[0])
        qual2.append(both_quals[1])
        title.append(read[0])

    # create the consensus reads
    # read1
    con_seq1, con_qual1 = combine_reads(seq1, qual1, min_qual)
    # read2
    con_seq2, con_qual2 = combine_reads(seq2, qual2, min_qual)

    # check the reads
    bad_seq1_count = consensus_checker(con_seq1, seq1, consen_cutoff)
    bad_seq2_count = consensus_checker(con_seq2, seq2, consen_cutoff)

    # generate the new name
    updated_name = title[0].split("\\")[0] + "con"
    name_record = {updated_name: title}

    # generate strings with both reads for output
    con_seq = con_seq1 + "_" + con_seq2
    con_qual = con_qual1 + "_" + con_qual2

    # create read tuple: (title, sequence, qual)
    updated_consol_read = (updated_name, con_seq, con_qual)

    return updated_consol_read, name_record, bad_seq1_count, bad_seq2_count
