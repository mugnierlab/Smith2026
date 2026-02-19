# Figure 2 F
# DNA double stranded breaks trigger mosaic VSG formation if homology is available
# Smith 2024

# read in data

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

# determine recombination length and plot

recomb_lenE9 <- all_E9 %>%
  mutate(recomb_len = target_end - target_start)

recomb_lenD3 <- all_D3_reseq %>%
  mutate(recomb_len = target_end - target_start)


all_recomb_sitesE9_D3 <- bind_rows(recomb_lenD3, recomb_lenE9)

all_recomb_sites_graph <- ggplot(all_recomb_sitesE9_D3) +
  geom_histogram(aes(x = recomb_len), binwidth = 5) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 20,
                                                     angle = 90,
                                                     hjust = 1,
                                                     vjust = 0.5),
                 axis.text.y = ggplot2::element_text(size = 20),
                 axis.title.y = ggplot2::element_text(size = 25),
                 axis.title.x = ggplot2::element_text(size = 25),
                 plot.background = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(),
                 panel.border = element_rect(fill = NA,colour = "grey50"),
                 axis.line = ggplot2::element_line(color = "black", size = 2),
                 axis.ticks.length = ggplot2::unit(0.25, "cm"),
                 legend.position = "right",
                 strip.text.x = ggplot2::element_text(size = 20),
                 legend.title = ggplot2::element_text(size = 15,
                                                      face = "bold",
                                                      hjust = 0.5),
                 legend.key = ggplot2::element_rect(fill = "white"),
                 legend.text = ggplot2::element_text(size = 15)) +
  ggplot2::labs(x = "Recomb Site Length", y = "Count") #+
coord_cartesian(xlim = c(0, 200))
