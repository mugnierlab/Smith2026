# Supplemental Figure 7A
# Smith 2025

# read in data
all_short_trunc <- get_all_CO_events("./ind_cross_over_events_truncations.csv", day_order = c("DOX", "RAD51KO"),
                                     genotype_order = c("AnTat.592.100c", "AnTat.592.200c", "AnTat.592.300c", "AnTat.592.400c", "AnTat.592.100b", "AnTat.592.200b", "AnTat.592.300b", "AnTat.592.400b", "NoGuide"), mouse_order = c("E4-C8", "Mun1"),
                                     primer_order = c("1F", "2F", "3F", "4F", "0R", "1R", "2R", "3R"))

all_short_trunc_filtered <- all_trunc %>%
  filter(mouse != "Mun1") %>%
  mutate(batch = case_when(
    str_detect(genotype, "c") ~ "bactch1",
    str_detect(genotype, "b") ~ "batch2",
    TRUE ~ "none"
  ))  %>%
  mutate(mouse = case_when(
    str_detect(genotype, "c") ~ "E4-C8a",
    str_detect(genotype, "b") ~ "E4-C8b",
    TRUE ~ mouse
  ))# %>%

all_short_trunc_filtered$genotype <- all_short_trunc_filtered$genotype %>%
  str_replace("c", "")

all_short_trunc_filtered$genotype <- all_short_trunc_filtered$genotype %>%
  str_replace("b", "")

all_short_trunc_filtered <- all_short_trunc_filtered %>%
  filter(type_of_event %in% c("target_to_donor", "donor_to_target", "donor_to_target_edge",
                              "target_to_donor_edge"))


all_long_trunc<- get_all_CO_events("./ind_cross_over_events_truncations_long.csv", day_order = c("DOX"),
                                   genotype_order = c("AnTat.592.Parental", "AnTat.592.rDNA", "AnTat.592.rDNA-noUTR-C", "AnTat.592.rDNA-500-C", "AnTat.592.rDNA-100-C"), mouse_order = c("E4-C8a", "E4-C8b"),
                                   primer_order = c("1F", "2F", "3F", "4F", "0R", "1R", "2R", "3R"))
all_long_trunc <- all_long_trunc %>%
  filter(type_of_event %in% c("target_to_donor", "donor_to_target", "donor_to_target_edge",
                              "target_to_donor_edge")) %>%
  mutate(batch = "batch0")

all_trunc_data <- bind_rows(all_long_trunc, all_short_trunc_filtered)

all_trunc_data <- qPCRr::reorder_samples(all_trunc_data, "genotype", c("AnTat.592.Parental", "AnTat.592.rDNA", "AnTat.592.rDNA-noUTR-C", "AnTat.592.rDNA-500-C", "AnTat.592.400", "AnTat.592.300", "AnTat.592.200", "AnTat.592.100", "AnTat.592.rDNA-100-C"))

# supplemental
graph_all_events(all_trunc_data, color_var = "donor_VSG", wrap_by = c("genotype", "mouse"), leg = "Donor VSG", adjust_xaxis = TRUE, multiple = FALSE, vert_num = 9, bin_width = 5) +
  ggplot2::geom_vline(aes(xintercept = c(86))) +
  ggplot2::geom_vline(aes(xintercept = c(1598))) +
  ggplot2::geom_vline(aes(xintercept = c(1079))) +
  ggplot2::geom_vline(aes(xintercept = c(693))) + 
  ggplot2::geom_vline(aes(xintercept = c(643))) +
  coord_cartesian(xlim = c(400, 1000))
