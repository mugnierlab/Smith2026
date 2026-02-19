# Smith2025
identifying mosaic VSGs\
please email jaclynsmith380@gmail.com and mmugnie1@jhu.edu with any questions

VSG-AMP-Seq: A detailed Readme can be found within the VSG-AMP-Seq Folder. A vignette is included.

Figures: all code used to generate figures is present in each subfolder. Relevant data needed to generate the figures is also present.

VSG_clustering: two alternative clustering methods were used to cluster VSGs into families. The code used to obtain these clusters and the resulting clustered outputs are in this folder.

VSgnome_analysis: Additional VSG-encoding genes were identified in the 2018 Lister Genome (Müller et al., 2018 Nature).
    ```
    genome_VSG_analysis.py
    ```
This script runs the pipeline - many large data files were organized into folders with paths relative to the working directory containing the scripts. Most data files are found on TriTrypDB.com. ABEX7_all_orfs.fa is from Beaver et al. 2024. All data files can be provided upon request, but are too large for upload onto github.

altVSG8_nanopore_analysis: scripts and data files to analyze VSG amplicons from muMT mice infected with Tb427VSG-8 expressing parasites. ```nanopore_mosaic_VSG_amplicon_pipleine.py``` 
identifies mosaic VSG reads from raw nanopore data against a specific donor VSG. 
```nanopore_mosaic_VSG_quant.py```
normalizes samples for read depth by counting all VSG-containing reads per sample.
Since mosaic nanopore reads were corrected by hand, all corrected reads are present here in 
```all_corrected_mos_ol.fasta```

aux_scripts: Additional scripts used to convert .ab1 files into fastqs and for FASTA processing

donor_intact_assay: Plasmidsaurus nanopore amplicon sequencing of the donor VSG amplicon, corresponding mosaic VSGs (.ab1 converted to fastq), and control mosaic VSGs from the parental line, Monomorphic Single Marker 427 1339 Cas9 TetR T7RNAP, to demonstrate VSG-228 is not present. (.ab1 converted to fastq & raw .ab1 files present).

genomic_sketch_R: the functions used to plot the genome in R. \
libraries required: tidyverse

nanopore_colony_consensus_builder: Plasmidsaurus fragmented raw nanopore VSG sequencing reads from mosaic VSG colonies (Original VSGs: Tb427VSG-8 and Tb1125VSG-73). 
```plasmid_nanopore_con_seq_pipeline.py```
determines the full length VSGs from a mixed colony with Plasmidsaurus nanopore reads using cd-hit. 
```package_list.txt``` specifies packages installed into the environment to run the script.

tissue_ORF_analysis: the python script used to quickly identify ORFs from Beaver et al. 2024 which were AnTat1.1 mosaics\
requires bowtie(1.2.3) \
The output of this analysis is formatted for R functions from VSG-AMP-Seq. Unknown mosaic VSGs with an identified portion matching AnTat1.1 were parsed by hand.

all_AnTat_mosaic_clones.fastq - a FASTQ file with all the full length VSG sequences from the individually isolated clones after guide induction in EATRO1125 cells targeting VSG AnTat1.1.

all_AnTat_mosaic_colonies_deathAssays.fastq - a FASTQ file with all the colony consensus sequences from mixed populations of mosaic parasites after guide induction in EATRO1125 cells targeting VSG AnTat1.1. Low base quality indicates a mixed population at that position.

all_Lister_cut_clones.fastq - a FASTQ file with all the full length VSG sequences from the individually isolated clones after guide induction in Lister427 cells targeting VSG-2.

all_VSG-73_mosaic_sequences.fasta - a FASTA file with all the full length consensus sequences & clonal sanger sequences from isolated mosaic VSGs. Further detials can be found in supplemental excel file 1

all_VSG-8_mosaic_sequences.fasta - a FASTA file with all the full length consensus sequences & clonal sanger sequences from isolated mosaic VSGs. Further detials can be found in supplemental excel file 1
