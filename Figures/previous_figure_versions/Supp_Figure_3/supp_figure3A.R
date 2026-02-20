# figure supplemental 3A
# Smith 2025

# read in data
all_D3_reseq <- get_all_CO_events("./ind_cross_over_events_D32.csv", day_order = c("DOX", "DMSO"),
                                  genotype_order = c("NoGuide", "Antat.141", "Antat.300R", "Antat.592",
                                                     "Antat.792", "Antat.909R", "Antat.1357"), mouse_order = c("D3"),
                                  primer_order = c("1F", "2F", "3F", "4F", "0R", "1R", "2R", "3R"))

all_D3_reseq <- all_D3_reseq %>%
  filter(type_of_event %in% c("target_to_donor", "donor_to_target", "target_to_donor_edge", "donor_to_target_edge"))

# summary of all recomb sites
all_events_graph_D3_2 <- graph_all_events(all_D3_reseq, color_var = "donor_VSG", wrap_by = "genotype", cut_line = TRUE,
                                          cut_line_names = c("NoGuide", "Antat.141", "Antat.300R", "Antat.592", "Antat.792",
                                                             "Antat.909R", "Antat.1357"), cut_line_locations = c(0, 242, 369, 693, 893, 979, 1458),
                                          cut_line_colors = c("black"), leg = "Donor VSG" , adjust_xaxis = TRUE, vert_num = 7, log_scale = FALSE, multiple = FALSE)+
  ggplot2::geom_vline(aes(xintercept = c(86))) +
  ggplot2::geom_vline(aes(xintercept = c(1598))) +
  ggplot2::geom_vline(aes(xintercept = c(1079)))

all_events_graph_D3_2