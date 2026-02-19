VSG-AMP-Seq Pipeline\
Author: Jaclyn Smith

Install Instructions\
Tested on: MacOS, MacOS with M1 Chip & Linux\
Install time: 15 minutes\
Create anaconda environment and install following packages:

```
  conda config --add subdirs osx-64 # needed for M1 chip
  conda install anaconda::python # 3.8.19
  conda install anaconda::biopython # 1.78
  conda config --add channels bioconda
  conda config --add channels conda-forge
  conda install bioconda::cutadapt (check this must be 3.2>) # 3.5
  conda install bioconda::trim-galore # v0.6.4
  conda install anaconda::pandas #1.2.1
  conda install conda-forge::python-levenshtein #0.25.1
  pip install progressbar # OR conda install conda-forge::progressbar #2.5
  conda install -c conda-forge regex #2.5.82
  conda install bioconda::bowtie #1.3.1
```

To run:\
  Required files:\
  \
    global_target.py - this file must be updated for the specific target in question:\
      full_length_target - a string with the full sequence from the splice leader to the 14-mer\
      alt_splicing - a dictionary, the observed alt splice events and/or errors detected near the\
          splice leader with full length sequences as the values. This is required for speeding up\
          the pipeline\
      protein_start - an int, the start position of the protein coding sequence. 0 based\
      antat_primers.txt - the name of the .txt file with primers and sequences\
    target.fasta - a FASTA with the full length sequence of the target from splice leader to 14mer\
    target_primers.txt - a tsv format file with primer#, primer sequence, and the predicted sequence of the anchored read\
    target_primers.fasta - a FASTA with primer names and the sequence names and the primers as the sequence\
              These names are used for the output files and should include the target_primer#\
    expt_barcodes.txt - a tsv format file with expt name - must be a 3 part name with experimental\
              conditions separated by underscores\
              examples Mouse_day_genotype, or cell-line_drug_target-cut-site\
    quant_R1s_bowtie.py - update the positions of the guides based within the full length target. 0 based\
  \
  Edit the main() within each module to update for the specific input files. Alternatively can edit the run_pipeline.py file from the vignette.

  The amount of time it takes to run these scripts depends on the data size, the number of mosaics within the sample, and the number of donor VSGs to check. As each of these increases, the length of the run also increases. To shorten the run time, the number of donor VSGs can be limited, reads can undergo consolidation, and samples can be parallelized. The pipeline has been optimized to run even large datasets off the NovaSeq6000 and should take less than a day to obtain demultiplexed files on a basic Mac laptop.
  I recommend using tmux or screen to set up a number of parallel instances of VSG_align. They can run without interfering with each other. If a file is running very slowly - there may be a way to filter out common unknown reads without testing against each potential donor VSG.

  genome_sort_by_primer:\
    genome_primer_sort(<VSG_library_of_interest.fa>, <target_primers.txt>)
  
\
  Index_addition:\
    index_add("<R1_file>", "<R2_file>",\
              "<I1_file>", "<I2_file>", "NovaSeq" or "MiSeq")\
\
  primer_sort:\
    sort_sequences("expt_code", "target_primers.fasta", "./index_reads/index_R1.fastq", "./index_reads/index_R2.fastq")\
\
  trim:
  * for this file you must remove &> if you are on linux and not mac!!!\
    primer_dict, prime_seq = aux_functions.read_in_primers("target_primers.txt")\
    trim("./sorted_reads/", primer_dict, <Target>) # Target is a String like "AnTat"\
    spacer_trim(primer_dict, <Target>)\
\
  demultiplex:\
    bc = barcodes("expt_barcodes.txt")\
    allowed_error = barcode_errors(bc)\
    demultiplex(bc, allowed_error)\
            if error because too many files are open, run in terminal:
            ```
            ulimit -n <num of files to open> 
            ```\
\
  consol_reads:\
    note files in the demulitplexed_reads folder must be gzipped!\
  consol_reads_updated: this is for larger files, needs use of a super computer and written for linux.
\
  vsg_align:\
    vsg_align(<source_folder>) # must have unzipped files; "demultiplexed_reads" or "consol_reads"\
\
  define_read_consensus:\
    no changes\
\
  identify_mosaics:\
    mosaic_find("target_primers.txt")\
\
  run mosaic_graphing_functions:\
  To run, must load libraries: devtools, tidyverse\
  run command:
  ```
install_github("JaclynSmith/qPCRr")
```
load library qPCRr \
  \
    get_all_CO_events(file_name, day_order, genotype_order, mouse_order, primer_order)\
      file_name = string - path to file\
      day_order = vector with names of second experimental attribute\
      genotype_order = vector with names of third experimental attribute\
      mouse_order = vector with names of first experimental attribute\
      primer_order = vector with names of the primers\
      # note if there are any missing samples, these need to also be absent from the vectors

  \
    graph_all_events(single_events) # this is the output of get_all_CO_events, many options are toggleable\
    All settings for each figure are in Figures section.

