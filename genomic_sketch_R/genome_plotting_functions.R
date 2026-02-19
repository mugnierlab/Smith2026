# functions for plotting genes and annotations to genome
# Smith 2024
#

add_genes <- function(gff_file, color, position_tsv, to_be_plotted = "Chr", plotting_ratio = conversion_factor,
                      separate_chrs = c("chr", "phase", "version"), sep = "_"){
  #read in positions of start of chr within plotting window
  positions <- readr::read_tsv(position_tsv, col_names = TRUE)
  gff <- readr::read_tsv(gff_file,
                         col_names = c("chr", "source", "type", "start", "end", "score", "strand",
                                       "phase_ignore", "attributes"))
  gff <- gff %>%
    tidyr::drop_na() %>%
    dplyr::filter(stringr::str_detect(chr, to_be_plotted))
  
  gff <- gff %>%
    dplyr::mutate(start_pos = start/plotting_ratio,
                  end_pos = end/plotting_ratio) %>%
    tidyr::separate(chr, separate_chrs, sep) %>%
    dplyr::select(-start, -end) %>%
    dplyr::mutate(multiplier = dplyr::case_when(
      strand == "-" ~ -1,
      strand == "+" ~ 1,
      TRUE ~ 0
    ))
  plotting_gff <- dplyr::left_join(gff, positions, by = c("chr", "phase"))
  
  for (i in 1:length(plotting_gff$chr)) {
    rect(xleft = plotting_gff$start_x[i] + plotting_gff$start_pos[i],
         xright = plotting_gff$start_x[i] + plotting_gff$end_pos[i],
         ybottom = plotting_gff$start_y[i],
         ytop = plotting_gff$start_y[i] + (0.10*plotting_gff$multiplier[i]),
         border = NA,
         col = color)
  }
}

add_annotations <- function(gff_file, color, position_tsv, symbol_print = "*", to_be_plotted = "Chr",
                            plotting_ratio = conversion_factor, separate_chrs = c("chr", "phase", "version"),
                            sep = "_", join_by = c("chr", "phase"), sep_bool = TRUE){
  #read in positions of start of chr within plotting window
  positions <- readr::read_tsv(position_tsv, col_names = TRUE)
  gff <- readr::read_tsv(gff_file,
                         col_names = c("chr", "source", "type", "start", "end", "score", "strand",
                                       "phase_ignore", "attributes"))
  # only plots genes off of a particular location
  gff <- gff %>%
    tidyr::drop_na() %>%
    dplyr::filter(stringr::str_detect(chr, to_be_plotted))
  
  if (sep_bool) {
    gff <- gff %>%
      dplyr::mutate(start_pos = start/plotting_ratio,
                    end_pos = end/plotting_ratio) %>%
      tidyr::separate(chr, separate_chrs, sep) %>%
      dplyr::select(-start, -end) %>%
      dplyr::mutate(multiplier = dplyr::case_when(
        strand == "-" ~ -1,
        strand == "+" ~ 1,
        strand == "." ~ 1,
        TRUE ~ 0
      ))
  }
  else {
    gff <- gff %>%
      dplyr::mutate(start_pos = start/plotting_ratio,
                    end_pos = end/plotting_ratio) %>%
      dplyr::select(-start, -end) %>%
      dplyr::mutate(multiplier = dplyr::case_when(
        strand == "-" ~ -1,
        strand == "+" ~ 1,
        strand == "." ~ 1,
        TRUE ~ 0
      ))
  }

  plotting_gff <- dplyr::left_join(gff, positions, by = join_by)

  for (i in 1:length(plotting_gff$chr)) {
    x_pos = mean(c(plotting_gff$start_x[i] + plotting_gff$start_pos[i],
                   plotting_gff$start_x[i] + plotting_gff$end_pos[i]))
    y_pos = plotting_gff$start_y[i] + (0.25*plotting_gff$multiplier[i])
    text(x_pos, y_pos, symbol_print, col = color, adj = 0.5)
  }
  
}
