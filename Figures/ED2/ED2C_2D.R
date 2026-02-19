# Extended Data Figure 2C & 2D
# Cas9 system and sgRNA design and analysis
# Smith 2026

# read in data

all_possible_recombs <- get_all_CO_events("./final_recomb_sites.csv", day_order = c("D1"),
                                          genotype_order = c("geno"),
                                          mouse_order = c("M1"),
                                          primer_order = c("primer"))

# ED 2C graph the recombination events
rep_donor_graph <- graph_all_events(all_possible_recombs, wrap_by = "genotype",
                                    color_var = "donor_VSG", leg = "Donor VSGs", multiple = TRUE,
                                    color_scheme = TRUE, adjust_xaxis = TRUE) +
  ggplot2::geom_vline(aes(xintercept = c(86))) +
  ggplot2::geom_vline(aes(xintercept = c(1598))) +
  ggplot2::geom_vline(aes(xintercept = c(1079)))


# ED 2D 
# determine recombination events per guide
rs_141 <- all_possible_recombs %>%
  filter(0 <= avg_pos & avg_pos <= 493) %>%
  group_by(donor_VSG) %>%
  summarise(count = n()) %>%
  add_row(donor_VSG = "VSG-4156", count = 0) %>%
  add_row(donor_VSG = "VSG-7358", count = 0) %>%
  mutate(guide = "G_243")


rs_300R <-all_possible_recombs %>%
  filter(369-250 <= avg_pos & avg_pos <= 369+250) %>%
  group_by(donor_VSG) %>%
  summarise(count = n()) %>%
  add_row(donor_VSG = "VSG-4156", count = 0) %>%
  add_row(donor_VSG = "VSG-7358", count = 0) %>%
  mutate(guide = "G_369R")

rs_592 <-all_possible_recombs %>%
  filter(694-250 <= avg_pos & avg_pos <= 694+250) %>%
  group_by(donor_VSG) %>%
  summarise(count = n()) %>%
  mutate(guide = "G_694")

rs_792 <-all_possible_recombs %>%
  filter(894-250 <= avg_pos & avg_pos <= 894+250) %>%
  group_by(donor_VSG) %>%
  summarise(count = n()) %>%
  mutate(guide = "G_894")

rs_909R <-all_possible_recombs %>%
  filter(978-250 <= avg_pos & avg_pos <= 978+250) %>%
  group_by(donor_VSG) %>%
  summarise(count = n()) %>%
  mutate(guide = "G_978R")

rs_1357 <-all_possible_recombs %>%
  filter(1459-250 <= avg_pos & avg_pos <= 1459+250) %>%
  group_by(donor_VSG) %>%
  summarise(count = n()) %>%
  mutate(guide = "G_1459")

recomb_stats_perGUide <- bind_rows(rs_141, rs_300R, rs_592, rs_792, rs_909R, rs_1357)

recomb_stats_perGUide <- qPCRr::reorder_samples(recomb_stats_perGUide, "guide", c("G_243",
                                                                                  "G_369R",
                                                                                  "G_694",
                                                                                  "G_894",
                                                                                  "G_978R",
                                                                                  "G_1459"))

recomb_guide_summary <- recomb_stats_perGUide %>%
  group_by(guide) %>%
  summarise(avg = mean(count), stdev = sd(count))


ggplot(recomb_stats_perGUide) +
  geom_jitter(aes(x = guide, y = count, color = donor_VSG), width = 0.25) +
  geom_crossbar(data = recomb_guide_summary, aes(x = guide, y = avg, ymin = avg, ymax = avg)) +
  geom_errorbar(data = recomb_guide_summary, aes(x = guide, y = avg, ymin = avg - stdev, ymax = avg + stdev)) +
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
  scale_color_manual(values = color_VSG_vector,
                     limits = c("VSG-228", "VSG-2986", "VSG-3110", "VSG-7358",
                                "VSG-4156", "VSG-1128", "VSG-479", "other VSGs"),
                     na.value = "grey80") #+
coord_cartesian(ylim = c(0,20))
