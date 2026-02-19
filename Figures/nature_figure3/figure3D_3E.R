# Figure 3D & 3E 
# Mosaic formation is dependent upon components of homologous recombination
# Smith 2026
#rad51 kos in vitro

all_invitroRad51 <- get_all_CO_events("./ind_cross_over_events_rad51_deepSeq.csv", day_order = c("WT", "RAD51KO", "B2KO"),
                                    genotype_order = c("AnTat.592.B6N1a", "AnTat.592.B6N1b", "AnTat.592.C4", "AnTat.592.J1339a", "AnTat.592.J1339b"), mouse_order = c("Mun1"),
                                    primer_order = c("1F", "2F", "3F", "4F", "0R", "1R", "2R", "3R"))
View(all_invitroRad51)

all_invitroRad51 <- all_invitroRad51 %>%
  filter(type_of_event %in% c("target_to_donor", "donor_to_target", "target_to_donor_edge", "donor_to_target_edge"))

# 3D graph
all_invitroRad51_graph <- graph_all_events(all_invitroRad51, color_var = "donor_VSG", wrap_by = c("day", "genotype"), cut_line = TRUE,
                                           cut_line_names = c("Antat.592"), cut_line_locations = c(693),
                                           cut_line_colors = c("black"),
                                         leg = "Donor VSG", adjust_xaxis = TRUE, multiple = FALSE, vert_num = 3, log_scale = FALSE) #+
  ggplot2::geom_vline(aes(xintercept = c(86))) +
  ggplot2::geom_vline(aes(xintercept = c(1598))) +
  ggplot2::geom_vline(aes(xintercept = c(1079)))

# 3E graph
summary_guides_rad51ko <- all_invitroRad51 %>%
  mutate(pos = case_when(
    genotype == "NoGuide" ~ 0,
    genotype == "Antat.141" ~ 242,
    genotype == "Antat.300R" ~ 369,
    genotype == "Antat.592" ~ 693,
    genotype == "Antat.792" ~ 893,
    genotype == "Antat.909R" ~ 979,
    genotype == "Antat.1357" ~ 1458,
    TRUE ~ 0
  )) %>%
  group_by(day, genotype) %>%
  filter((avg_pos >= pos - 250) | (avg_pos <= pos + 250)) %>%
  summarise(count = n()) %>%
  filter(genotype != "NoGuide") %>%
  ungroup() %>%
  add_row(day = "RAD51", genotype = "AnTat.592.B1N2a", count = 0) %>%
  add_row(day = "RAD51", genotype = "AnTat.592.B1N2b", count = 0) %>%
  add_row(day = "RAD51", genotype = "AnTat.592.C1", count = 0) %>%
  add_row(day = "RAD51", genotype = "AnTat.592.C2", count = 0) %>%
  add_row(day = "RAD51", genotype = "AnTat.592.C4", count = 0)
  


rad51_R1 <- read_tsv("output_R1s_rad51_invitro.tsv", col_names = c("file_info", "count"))

rad51_R1$file_info <- str_replace(rad51_R1$file_info,"./sam_alignments/", "")
rad51_R1_sum <- rad51_R1 %>%  
  separate(file_info, sep = "_", into = c("mouse", "day", "genotype", "barcode", "primer", "extra" )) %>%
  select(-extra) %>%
  group_by(genotype) %>%
  summarise(total_reads = sum(count))

rad51_final <- full_join(summary_guides_rad51ko, rad51_R1_sum, by = "genotype")

rad51_final_sum <- rad51_final %>%
  mutate(ratio = total_reads/6014614, expt = "rad51KO") %>%
  mutate(multiplier_account = count/ratio)

rad51_final_sum$genotype = str_replace(rad51_final_sum$genotype, "AnTat.", "")
rad51_final_sum <- rad51_final_sum %>%
  separate_wider_delim(genotype, names=c("genotype", "clone"), delim = ".")


rad51_ratios <- rad51_final_sum %>%
  mutate(genotype = case_when(
    genotype == "141" ~ "guide_243",
    genotype == "300R" ~ "guide_369*",
    genotype == "592" ~ "guide_694",
    genotype == "792" ~ "guide_894",
    genotype == "909R" ~ "guide_978*",
    genotype == "1357" ~ "guide_1459",
    TRUE ~ ""
  ))

rad51_ratios <- rad51_ratios %>%
  qPCRr::reorder_samples("day", c("WT", "RAD51", "B2KO"))

sum_rad51_ratios <- rad51_ratios %>%
  group_by(day) %>%
  summarise(avg = mean(multiplier_account), stdev = sd(multiplier_account))

# res_aov <- aov(multiplier_account ~ day, data = rad51_ratios)
# summary(res_aov)
# TukeyHSD(res_aov)


ggplot(rad51_ratios) +
  geom_point(aes(x = day, y = multiplier_account, color = clone), position = position_dodge(width = 0.5)) +
  geom_crossbar(data = sum_rad51_ratios, aes(x = day, y = avg, ymin = avg, ymax = avg)) +
  geom_errorbar(data = sum_rad51_ratios, aes(x = day, y = avg, ymin = avg - stdev, ymax = avg + stdev)) +
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
  ggplot2::labs(x = "Guides", y = "normalized mosaic count", color = "Clone")
