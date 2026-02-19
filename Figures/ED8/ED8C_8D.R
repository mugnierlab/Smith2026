# Extended Data Figure 8C&8D
# Rad51 KOs in mumT mice 
# Smith 2026

# read in data
all_mouseRad51 <- get_all_CO_events("./ind_cross_over_events_muMTRad51.csv", day_order = c("D7", "D8", "D14"),
                            genotype_order = c("muMT"), mouse_order = c("MR2", "MR3", "MR4", "MR5"),
                            primer_order = c("1F", "2F", "3F", "4F", "0R", "1R", "2R", "3R"))
View(all_mouseRad51)

all_mouseRad51 <- all_mouseRad51 %>%
  filter(type_of_event %in% c("target_to_donor", "donor_to_target", "target_to_donor_edge", "donor_to_target_edge"))

# graph Figure 8C
all_mouseRad51_graph <- graph_all_events(all_mouseRad51, color_var = "donor_VSG", wrap_by = "mouse", 
                                        leg = "Donor VSG", adjust_xaxis = TRUE, multiple = TRUE, vert_num = 7, log_scale = FALSE) #+
  ggplot2::geom_vline(aes(xintercept = c(86))) +
  ggplot2::geom_vline(aes(xintercept = c(1598))) +
  ggplot2::geom_vline(aes(xintercept = c(1079)))

# read in data
rad51_samples_consol_totals <- readr::read_csv("./target_alignment_stats_muMTRad51.csv")

reads_per_sample_rad51 <- rad51_samples_consol_totals %>%
  group_by(mouse, day, genotype) %>%
  summarise(total_reads = sum(target_R1))

# norm to MR5; D8
reads_per_sample_rad51 <- reads_per_sample_rad51 %>%
  mutate(normalize = 751068) %>%
  mutate(ratio = total_reads/normalize)


rad51_quant <- all_mouseRad51 %>%
  ungroup() %>%
  group_by(mouse, day, genotype) %>%
  summarise(total = n()) %>%
  ungroup() %>%
  dplyr::add_row(mouse = "MR3", day = "D7", genotype = "muMT", total = 0)


 rad51_summary_stats <- rad51_quant %>%
  full_join(reads_per_sample_rad51, by = c("mouse", "day", "genotype")) %>%
  mutate(norm_total = total/ratio)

rad51_summary_stats <- rad51_summary_stats %>%
  mutate(day_groups = case_when(
    day == "D7" ~ "early",
    day == "D8" ~ "early",
    day == "D14" ~ "late"
  ))


rad_51_summary <- rad51_summary_stats %>%
  ungroup() %>%
  group_by(day_groups) %>%
  summarise(avg = mean(norm_total), stdev = sd(norm_total))
  


# graph plot 8D
ggplot(data = rad51_summary_stats) +
  geom_jitter(aes(x = day_groups, y = norm_total)) +
  geom_crossbar(data = rad_51_summary, aes(x = day_groups, y = avg, ymin = avg, ymax = avg)) +
  geom_errorbar(data = rad_51_summary, aes(x = day_groups, y = avg, ymin = avg - stdev, ymax = avg + stdev)) +
  coord_cartesian(ylim = c(0, NA)) + 
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
  ggplot2::labs(x = "Day", y = "Count", fill = "Experiment")

stats::shapiro.test(rad51_summary_stats$norm_total)

rad51_sum_early <- rad51_summary_stats %>%
  filter(day_groups == "early")
rad51_sum_late <- rad51_summary_stats %>%
  filter(day_groups == "late")

wilcox.test(rad51_sum_early$norm_total, rad51_sum_late$norm_total, paired=TRUE)
