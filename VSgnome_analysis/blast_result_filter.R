# filter blast results to find unknown proteins which can reasonably be called VSGs

VSG_blast_results <- readr::read_tsv("../VSGnome_analysis_jaclyn/unknown_protein_VSGBlast/blastResults.tsv",
                                     col_names = c("unknown", "putative_VSG", "percent_identity", "number_identical",
                                                   "length", "mismatch", "gap_open", "unknown_start", "unknown_end",
                                                   "unknown_length", "VSG_start", "VSG_end", "VSG_length", "evalue",
                                                   "bitscore"))
VSG_blast_trim <- VSG_blast_results %>%
  # filter out those where the full unknown length isn't covered by the query
  filter(length/unknown_length*100 >= 80) %>%
  # filter out those with a bitscore on lower end - some only partially contain VSGs
  filter(bitscore > 500)

# 2187 recovered VSG sequences
# 42 VSGs from Trinity assemblies from VSG-seq
# 


readr::write_tsv(VSG_blast_trim, "../VSGnome_analysis_jaclyn/unknown_protein_VSGBlast/blastResults_filtered.tsv")
