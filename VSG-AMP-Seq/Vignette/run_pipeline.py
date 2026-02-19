import genome_sort_by_primer
import index_addition
import primer_sort
import trim
import demultiplex
import vsg_align
import define_read_consensus
import identify_mosaics
import aux_functions
import os


def main():
    # run the functions
    genome_sort_by_primer.genome_primer_sort("VSGnome.fa", "antat_primers.txt")
    index_addition.index_add("test_R1.fastq", "test_R2.fastq",
                             "test_I1.fastq", "test_I2.fastq", "MiSeq")
    primer_sort.sort_sequences("TEST", "antat_primers.fasta", "./index_reads/index_R1.fastq", "./index_reads/index_R2.fastq")
    primer_dict, prime_seq = aux_functions.read_in_primers("antat_primers.txt")
    trim.trim("./sorted_reads/", primer_dict, "AnTat")
    trim.spacer_trim(primer_dict, "AnTat")
    bc = demultiplex.barcodes("test_barcodes.txt")
    allowed_error = demultiplex.barcode_errors(bc)
    demultiplex.demultiplex(bc, allowed_error)
    os.system("mkdir ./demultiplexed_reads/unknown_files")
    os.system("mv ./demultiplexed_reads/unassign* demultiplexed_reads/unknown_files")
    vsg_align.vsg_align("demultiplexed_reads")
    define_read_consensus.consensus_builder()
    identify_mosaics.mosaic_find("antat_primers.txt", "VSGnome.fa")




if __name__ == '__main__':
    main()
