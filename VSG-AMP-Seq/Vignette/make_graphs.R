indv <- get_all_CO_events("./ind_cross_over_events.csv", day_order = c("DOX"),
                          genotype_order = c("AnTat.592", "AnTat.592.BES7", "AnTat.592.Minichr", "AnTat.592.Parental", "AnTat.592.rDNA", "AnTat.592.Tubulin"),
                          mouse_order = c("E9", "D3", "D2-D3", "C7-B7"),
                          primer_order = c("1F", "2F", "3F", "4F", "0R", "1R", "2R", "3R"))

indv <- indv %>%
  filter(type_of_event %in% c("target_to_donor", "donor_to_target", "donor_to_target_edge",
                              "target_to_donor_edge"))

graph_all_events(indv, color_var = "donor_VSG", wrap_by = "genotype", leg = "Donor VSG", adjust_xaxis = TRUE, multiple = FALSE) +
  ggplot2::geom_vline(aes(xintercept = c(86))) +
  ggplot2::geom_vline(aes(xintercept = c(1598))) +
  ggplot2::geom_vline(aes(xintercept = c(1079))) +
  ggplot2::geom_vline(aes(xintercept = c(693)))
