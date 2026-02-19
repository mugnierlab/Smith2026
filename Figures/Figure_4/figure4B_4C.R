# Figure 4 B + C
# Donor VSGs are accessible throughout the genome
# Smith 2024

# read in data
all_kev <- get_all_CO_events("./ind_cross_over_events_kev.csv", day_order = c("DOX"),
                             genotype_order = c("NoGuide.Parental", "AnTat.592.Parental", "AnTat.592.Tubulin", "AnTat.592.rDNA", "AnTat.592.Minichr", "AnTat.592.BES7"), mouse_order = c("D2-1E9", "D2-D3", "C7-B7"),
                             primer_order = c("1F", "2F", "3F", "4F", "0R", "1R", "2R", "3R"))


all_kev_batchC7B7 <- all_kev %>%
  filter(type_of_event %in% c("target_to_donor", "donor_to_target", "donor_to_target_edge",
                              "target_to_donor_edge")) %>%
  filter(mouse %in% c("C7-B7", "D2-1E9"))

all_kev_batchD2D3 <- all_kev %>%
  filter(type_of_event %in% c("target_to_donor", "donor_to_target", "donor_to_target_edge",
                              "target_to_donor_edge")) %>%
  filter(mouse %in% c("D2-D3"))


# Figure 4B

graph_all_events(all_kev_batchD2D3, color_var = "donor_VSG", wrap_by = "genotype", leg = "Donor VSG", adjust_xaxis = TRUE, multiple = FALSE) +
  ggplot2::geom_vline(aes(xintercept = c(86))) +
  ggplot2::geom_vline(aes(xintercept = c(1598))) +
  ggplot2::geom_vline(aes(xintercept = c(1079))) +
  ggplot2::geom_vline(aes(xintercept = c(693)))

# Figure 4C

all_kev$mouse <- as.character(all_kev$mouse)
kev_summary <- all_kev %>%
  mutate(mouse = case_when(
    mouse == "D2-1E9" ~ "C7-B7",
    TRUE ~ mouse)) %>%
  group_by(mouse, day, genotype, donor_VSG) %>%
  filter((avg_pos >= 693 - 250) & (avg_pos <= 693 + 250)) %>%
  filter(type_of_event %in% c("target_to_donor", "donor_to_target", "target_to_donor_edge", "donor_to_target_edge")) %>%
  summarise(count = n()) %>%
  filter(genotype != "NoGuide.Parental") %>%
  filter(donor_VSG %in% c("VSG-2986", "VSG-228"))

kev_stats <- read_tsv("output_kev.tsv", col_names = c("file_info", "count"))
kev_stats$file_info <- str_replace(kev_stats$file_info,"./sam_alignments/", "")
kev_sum <- kev_stats %>%  
  separate(file_info, sep = "_", into = c("mouse", "day", "genotype", "barcode", "primer", "extra" )) %>%
  select(-extra) %>%
  group_by(mouse, genotype) %>%
  summarise(total_reads = sum(count)) %>%
  filter(genotype != "Antat.592") %>%
  mutate(mouse = case_when(
    mouse == "D2-1E9" ~ "C7-B7",
    TRUE ~ mouse
  ))


kev_quant_final <- full_join(kev_summary, kev_sum, by = c("mouse", "genotype")) %>%
  mutate(ratio = count/(total_reads/1407478)) %>%
  ungroup %>%
  add_row(mouse = "C7-B7", day = "DOX", genotype = "AnTat.592.Parental", donor_VSG = "VSG-228", count = 0, total_reads = 1, ratio = 0)
kev_quant_final <- qPCRr::reorder_samples(kev_quant_final, "genotype", c("AnTat.592.Parental", "AnTat.592.BES7", "AnTat.592.Minichr",
                                                                         "AnTat.592.rDNA", "AnTat.592.Tubulin"))
kev_quant_final <- qPCRr::reorder_samples(kev_quant_final, "donor_VSG", c("VSG-2986", "VSG-228"))

kev_ratios_sum <- kev_quant_final %>%
  ungroup() %>%
  group_by(genotype, donor_VSG) %>%
  summarise(avg = mean(ratio), stdev = sd(ratio))

write_csv(kev_quant_final, "kevin_quant.csv")

ggplot(kev_quant_final) +
  #ggplot() +
  geom_point(aes(x = genotype, y = ratio, color = mouse), position = position_dodge(width = 0.5)) +
  facet_wrap(~donor_VSG) +
  geom_crossbar(data = kev_ratios_sum, aes(x = genotype, y = avg, ymin = avg, ymax = avg)) +
  geom_errorbar(data = kev_ratios_sum, aes(x = genotype, y = avg, ymin = avg - stdev, ymax = avg + stdev)) +
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
  ggplot2::labs(x = "Insert Location", y = "log(normalized mosaic count)", color = "Clone") #+
# ggplot2::coord_cartesian(ylim = c(0, 4.5))

res_aov <- aov(ratio ~ interaction(genotype, donor_VSG), data = kev_quant_final)
summary(res_aov)
TukeyHSD(res_aov)
