# Extended Data Figure 4B
# DNA breaks result in switchers if a homologous donor is not present
# Smith 2026

# Clone 1 is D3, Clone 2 is E9

#read in data:
all_E9 <- get_all_CO_events("./ind_cross_over_events_E9.csv", day_order = c("DOX", "DMSO"),
                            genotype_order = c("NoGuide", "Antat.141", "Antat.300R", "Antat.592",
                                               "Antat.792", "Antat.909R", "Antat.1357"), mouse_order = c("E9"),
                            primer_order = c("1F", "2F", "3F", "4F", "0R", "1R", "2R", "3R"))

all_D3_reseq <- get_all_CO_events("./ind_cross_over_events_D32.csv", day_order = c("DOX", "DMSO"),
                                  genotype_order = c("NoGuide", "Antat.141", "Antat.300R", "Antat.592",
                                                     "Antat.792", "Antat.909R", "Antat.1357"), mouse_order = c("D3"),
                                  primer_order = c("1F", "2F", "3F", "4F", "0R", "1R", "2R", "3R"))

all_E9 <- all_E9 %>%
  filter(type_of_event %in% c("target_to_donor", "donor_to_target", "target_to_donor_edge", "donor_to_target_edge")) %>%
  filter(day != "DMSO")

all_D3_reseq <- all_D3_reseq %>%
  filter(type_of_event %in% c("target_to_donor", "donor_to_target", "target_to_donor_edge", "donor_to_target_edge")) %>%
  filter(day != "DMSO")

# graph rep sequences
all_events_graph_D3_2 <- graph_all_events(all_D3_reseq, color_var = "donor_VSG", wrap_by = "genotype", cut_line = TRUE,
                                          cut_line_names = c("NoGuide", "Antat.141", "Antat.300R", "Antat.592", "Antat.792",
                                                             "Antat.909R", "Antat.1357"), cut_line_locations = c(0, 242, 369, 693, 893, 979, 1458),
                                          cut_line_colors = c("black"), leg = "Donor VSG" , adjust_xaxis = TRUE, vert_num = 7, log_scale = FALSE, multiple = FALSE,
                                          summary = TRUE)

all_events_graph_E9 <- graph_all_events(all_E9, color_var = "donor_VSG", wrap_by = "genotype", cut_line = TRUE,
                                          cut_line_names = c("NoGuide", "Antat.141", "Antat.300R", "Antat.592", "Antat.792",
                                                             "Antat.909R", "Antat.1357"), cut_line_locations = c(0, 242, 369, 693, 893, 979, 1458),
                                          cut_line_colors = c("black"), leg = "Donor VSG" , adjust_xaxis = TRUE, vert_num = 7, log_scale = FALSE, multiple = FALSE,
                                          summary = TRUE)
