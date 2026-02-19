# Supplemental Figure 6
# The Lister 427 VSG Annotated Genome
# Smith 2025

plot.new(); plot.window(xlim = c(0, 8.5), ylim = c(0, 11))

# plot the genome scaffold

# obtains coordinates for each chromosome. need the full length of the things to plot. can be 
# obtained any way
coords <- readr::read_tsv("../VSGnome_analysis_jaclyn/genome_2018_coords.tsv", col_names = c("chr", "length"))

# selects only the core chromosomes for plotting, 1:700000 is selected as a conversion factor
# the length of the core chrs is determined here
conversion_factor = 700000
core_coords <- coords %>%
  filter(str_detect(chr, "core")) %>%
  mutate(norm_len = length/conversion_factor) %>%
  separate(chr, c("chr", "phase", "version"), "_")

# to deal with phasing and offsetting chrs such that they are drawn correctly
# which ever phase is longer is selected for the front length
front_both <- coords %>%
  filter(str_detect(chr, "_5")) %>%
  mutate(norm_len = length/conversion_factor) %>%
  filter(!str_detect(chr, "unitig")) %>%
  separate(chr, c("chr", "phase", "version"), "_") %>%
  group_by(chr) %>%
  slice(which.max(norm_len)) %>%
  select(chr, norm_len) %>%
  rename("front_norm_len" = "norm_len")

# determine the start positions based on the shift
core_coords_adjusted_start <- core_coords %>%
  left_join(front_both, by = "chr") %>%
  mutate(front_norm_len = replace_na(front_norm_len, 0)) %>%
  mutate(adjust_len = norm_len + front_norm_len)

# calculate the phase lengths
front_A <- coords %>%
  filter(str_detect(chr, "5A")) %>%
  mutate(norm_len = length/conversion_factor) %>%
  separate(chr, c("chr", "phase", "version"), "_")

front_B <-coords %>%
  filter(str_detect(chr, "5B")) %>%
  mutate(norm_len = length/conversion_factor) %>%
  separate(chr, c("chr", "phase", "version"), "_")

back_A <- coords %>%
  filter(str_detect(chr, "3A")) %>%
  mutate(norm_len = length/conversion_factor) %>%
  separate(chr, c("chr", "phase", "version"), "_")

back_B <- coords %>%
  filter(str_detect(chr, "3B")) %>%
  mutate(norm_len = length/conversion_factor) %>%
  separate(chr, c("chr", "phase", "version"), "_")

# store the positions of the chrs for later gene plotting
coordinate_starts <- tibble(
  chr = character(),
  phase = character(),
  start_x = numeric(),
  start_y = numeric()
)

# plot the chromosomes
for (i in 1:length(core_coords_adjusted_start$chr)) {
  
  print(core_coords_adjusted_start$chr[i])
  
  # print the core regions with offset for the longest 5' phased region
  lines(c(1 + core_coords_adjusted_start[i,]$front_norm_len, 1 + core_coords_adjusted_start[i,]$adjust_len), c(i, i))
  text(0.5, i, core_coords_adjusted_start$chr[i], srt = 90, cex = 0.5)
  coordinate_starts <- coordinate_starts %>%
    add_row(chr = core_coords_adjusted_start$chr[i],
            phase = "core",
            start_x = 1 + core_coords_adjusted_start[i,]$front_norm_len,
            start_y = i)
  
  # print the 5' scaffold and the phase for 5A
  if (core_coords_adjusted_start$chr[i] %in% front_A$chr) {
    front_A_length <- front_A %>%
      filter(chr == core_coords_adjusted_start$chr[i])
    
    lines(c(1 + core_coords_adjusted_start$front_norm_len[i] - front_A_length$norm_len,
            1 + core_coords_adjusted_start$front_norm_len[i]), c(i + 0.25, i + 0.25))
    lines(c(1 + core_coords_adjusted_start$front_norm_len[i], 1 + core_coords_adjusted_start$front_norm_len[i]), c(i, i+0.25))
    coordinate_starts <- coordinate_starts %>%
      add_row(chr = core_coords_adjusted_start$chr[i],
              phase = "5A",
              start_x = 1 + core_coords_adjusted_start$front_norm_len[i] - front_A_length$norm_len,
              start_y = i + 0.25)
  }
  
  # print the 5' scaffold and the phase 5B
  if (core_coords_adjusted_start$chr[i] %in% front_B$chr) {
    front_B_length <- front_B %>%
      filter(chr == core_coords_adjusted_start$chr[i])
    
    lines(c(1 + core_coords_adjusted_start$front_norm_len[i] - front_B_length$norm_len,
            1 + core_coords_adjusted_start$front_norm_len[i]), c(i - 0.25, i - 0.25))
    lines(c(1 + core_coords_adjusted_start$front_norm_len[i], 1 + core_coords_adjusted_start$front_norm_len[i]), c(i, i - 0.25))
    coordinate_starts <- coordinate_starts %>%
      add_row(chr = core_coords_adjusted_start$chr[i],
              phase = "5B",
              start_x = 1 + core_coords_adjusted_start$front_norm_len[i] - front_B_length$norm_len,
              start_y = i - 0.25)
  }
  
  # print the 3' scaffold and the phase 3A
  if (core_coords_adjusted_start$chr[i] %in% back_A$chr) {
    back_A_length <- back_A %>%
      filter(chr == core_coords_adjusted_start$chr[i])
    
    lines(c(1 + core_coords_adjusted_start$front_norm_len[i] + core_coords_adjusted_start$norm_len[i],
            1 + core_coords_adjusted_start$front_norm_len[i] + core_coords_adjusted_start$norm_len[i] + back_A_length$norm_len),
          c(i + 0.25, i + 0.25))
    lines(c(1 + core_coords_adjusted_start$front_norm_len[i] + core_coords_adjusted_start$norm_len[i],
            1 + core_coords_adjusted_start$front_norm_len[i] + core_coords_adjusted_start$norm_len[i]),
          c(i, i + 0.25))
    coordinate_starts <- coordinate_starts %>%
      add_row(chr = core_coords_adjusted_start$chr[i],
              phase = "3A",
              start_x = 1 + core_coords_adjusted_start$front_norm_len[i] + core_coords_adjusted_start$norm_len[i],
              start_y = i + 0.25)
  }
  
  # print the 3' scaffold and the phase 3B
  if (core_coords_adjusted_start$chr[i] %in% back_B$chr) {
    back_B_length <- back_B %>%
      filter(chr == core_coords_adjusted_start$chr[i])
    
    lines(c(1 + core_coords_adjusted_start$front_norm_len[i] + core_coords_adjusted_start$norm_len[i],
            1 + core_coords_adjusted_start$front_norm_len[i] + core_coords_adjusted_start$norm_len[i] + back_B_length$norm_len),
          c(i - 0.25, i - 0.25))
    lines(c(1 + core_coords_adjusted_start$front_norm_len[i] + core_coords_adjusted_start$norm_len[i],
            1 + core_coords_adjusted_start$front_norm_len[i] + core_coords_adjusted_start$norm_len[i]),
          c(i, i - 0.25))
    coordinate_starts <- coordinate_starts %>%
      add_row(chr = core_coords_adjusted_start$chr[i],
              phase = "3B",
              start_x = 1 + core_coords_adjusted_start$front_norm_len[i] + core_coords_adjusted_start$norm_len[i],
              start_y = i - 0.25)
  }
}


add_genes("../VSGnome_analysis_jaclyn/genome_2018_associatedFiles/TriTrypDB-66_TbruceiLister427_2018_VSGs.gff",
          "gray40",
          "../VSGnome_analysis_jaclyn/chr_coord_starts.tsv")
add_genes("../VSGnome_analysis_jaclyn/TriTrypDB-66_TbruceiLister427_2018_unknown_filtered.gff",
          "sandybrown",
          "../VSGnome_analysis_jaclyn/chr_coord_starts.tsv")

add_genes("../VSGnome_analysis_jaclyn/vsg_2986_family.gff",
          "purple",
          "../VSGnome_analysis_jaclyn/chr_coord_starts.tsv")

add_annotations("../VSGnome_analysis_jaclyn/vsg_2986_family.gff",
                "purple",
                "../VSGnome_analysis_jaclyn/chr_coord_starts.tsv", symbol_print = "*")

add_annotations("../VSGnome_analysis_jaclyn/insertion_site_gff/rDNA_spacer_inserts.gff",
                "blue",
                "../VSGnome_analysis_jaclyn/chr_coord_starts.tsv", symbol_print = "^")
add_annotations("../VSGnome_analysis_jaclyn/insertion_site_gff/tubulin_spacer_inserts.gff",
                "green",
                "../VSGnome_analysis_jaclyn/chr_coord_starts.tsv", symbol_print = "^")
# add scale bar
lines(c(1 + 5.5,
        1 + 5.5 + 500000/conversion_factor),
      c(1, 1))
text(6.85, 1.3, "0.5Mb", cex = 0.5)

# unitigs plotting

plot.new(); plot.window(xlim = c(0, 8.5), ylim = c(0, 11))

# store the positions of the chrs for later gene plotting
coordinate_starts <- tibble(
  chr = character(),
  phase = character(),
  start_x = numeric(),
  start_y = numeric()
)

unitig_list = c("unitig_2196_Tb427v10", "unitig_2457_Tb427v10", "unitig_1884_Tb427v10",
                "unitig_2456_Tb427v10", "unitig_1879_Tb427v10",
                "unitig_2198_Tb427v10", "unitig_1878_Tb427v10", "unitig_1881_Tb427v10",
                "unitig_223_Tb427v10", "unitig_252_Tb427v10")

unitig_plots <- coords %>%
  filter(chr %in% unitig_list) %>%
  #mutate(norm_len = length/conversion_factor)
  mutate(norm_len = length/10000)


for (i in 1:length(unitig_plots$chr)) {
  
  print(unitig_plots$chr[i])
  
  # print the core regions with offset for the longest 5' phased region
  lines(c(1.5, 1.5 + unitig_plots[i,]$norm_len), c(i, i))
  name = str_replace(unitig_plots$chr[i], "_Tb427v10", "")
  text(0.5, i, name, cex = 0.5)
  coordinate_starts <- coordinate_starts %>%
    add_row(chr = unitig_plots$chr[i],
            phase = "unitig",
            start_x = 1.5,
            start_y = i)
}

add_genes("../VSGnome_analysis_jaclyn/genome_2018_associatedFiles/TriTrypDB-66_TbruceiLister427_2018_VSGs.gff",
          "gray40",
          "../VSGnome_analysis_jaclyn/unitigs.tsv", plotting_ratio = 10000,
          to_be_plotted = "^unitig", join_by = c("chr"), sep_bool = FALSE)
add_genes("../VSGnome_analysis_jaclyn/TriTrypDB-66_TbruceiLister427_2018_unknown_filtered.gff",
          "sandybrown",
          "../VSGnome_analysis_jaclyn/unitigs.tsv", plotting_ratio = 10000,
          to_be_plotted = "^unitig", join_by = c("chr"), sep_bool = FALSE)
add_genes("../VSGnome_analysis_jaclyn/vsg_2986_family.gff",
          "purple",
          "../VSGnome_analysis_jaclyn/unitigs.tsv", plotting_ratio = 10000,
          to_be_plotted = "^unitig", join_by = c("chr"), sep_bool = FALSE)

add_annotations("../VSGnome_analysis_jaclyn/vsg_2986_family.gff",
                "purple",
                "../VSGnome_analysis_jaclyn/unitigs.tsv", symbol_print = "*", sep_bool = FALSE,
                join_by = c("chr"), to_be_plotted = "^unitig", plotting_ratio = 10000)

add_annotations("../VSGnome_analysis_jaclyn/insertion_site_gff/rDNA_spacer_inserts.gff",
                "blue",
                "../VSGnome_analysis_jaclyn/unitigs.tsv", symbol_print = "^", sep_bool = FALSE,
                join_by = c("chr"), to_be_plotted = "^unitig", plotting_ratio = 10000)
# add scale bar
lines(c(1 + 0,
        1 + 0 + 50000/10000),
      c(3, 3))
text(4, 1.3, "50kb", cex = 0.5)
