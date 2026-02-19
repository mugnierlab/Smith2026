# Figure 2 A + B
# DNA double stranded breaks trigger mosaic VSG formation if homology is available
# Smith 2026

# Panel A representative sequence is partial E9 clone, but D3 code is here also:
# These are also Extended Data Figure 3

#read in data:
all_E9 <- get_all_CO_events("./ind_cross_over_events_E9.csv", day_order = c("DOX", "DMSO"),
                            genotype_order = c("NoGuide", "Antat.141", "Antat.300R", "Antat.592",
                                               "Antat.792", "Antat.909R", "Antat.1357"), mouse_order = c("E9"),
                            primer_order = c("1F", "2F", "3F", "4F", "0R", "1R", "2R", "3R"))

all_D3_reseq <- get_all_CO_events("./ind_cross_over_events_D32.csv", day_order = c("DOX", "DMSO"),
                                  genotype_order = c("NoGuide", "Antat.141", "Antat.300R", "Antat.592",
                                                     "Antat.792", "Antat.909R", "Antat.1357"), mouse_order = c("D3"),
                                  primer_order = c("1F", "2F", "3F", "4F", "0R", "1R", "2R", "3R"))

# graph events
all_E9 <- all_E9 %>%
  filter(type_of_event %in% c("target_to_donor", "donor_to_target", "target_to_donor_edge", "donor_to_target_edge")) %>%
  filter(day != "DMSO")

# summary by guide position
all_events_graph_E9 <- graph_all_events(all_E9, color_var = "donor_VSG", wrap_by = "genotype", cut_line = TRUE,
                                             cut_line_names = c("NoGuide", "Antat.141", "Antat.300R", "Antat.592", "Antat.792",
                                                                "Antat.909R", "Antat.1357"), cut_line_locations = c(0, 242, 369, 693, 893, 979, 1458),
                                             cut_line_colors = c("black"), leg = "Donor VSG", adjust_xaxis = TRUE, multiple = FALSE, vert_num = 7, log_scale = FALSE) # +
  ggplot2::geom_vline(aes(xintercept = c(86))) +
  ggplot2::geom_vline(aes(xintercept = c(1598))) +
  ggplot2::geom_vline(aes(xintercept = c(1079)))

all_events_graph_E9

# summary of all recomb sites 
all_events_graph_E9_summary <- graph_all_events(all_E9, color_var = "donor_VSG", wrap_by = "type_of_event", cut_line = TRUE,
                                        cut_line_names = c("NoGuide", "Antat.141", "Antat.300R", "Antat.592", "Antat.792",
                                                           "Antat.909R", "Antat.1357"), cut_line_locations = c(0, 242, 369, 693, 893, 979, 1458),
                                        cut_line_colors = c("black"), leg = "Donor VSG", adjust_xaxis = TRUE, multiple = FALSE, vert_num = 7, log_scale = FALSE,
                                        summary = TRUE) +
  ggplot2::geom_vline(aes(xintercept = c(86))) +
  ggplot2::geom_vline(aes(xintercept = c(1598))) +
  ggplot2::geom_vline(aes(xintercept = c(1079)))


all_events_graph_E9_summary

all_D3_reseq <- all_D3_reseq %>%
  filter(type_of_event %in% c("target_to_donor", "donor_to_target", "target_to_donor_edge", "donor_to_target_edge"))

# summary of all recomb sites
all_events_graph_D3_2 <- graph_all_events(all_D3_reseq, color_var = "donor_VSG", wrap_by = "genotype", cut_line = TRUE,
                                          cut_line_names = c("NoGuide", "Antat.141", "Antat.300R", "Antat.592", "Antat.792",
                                                             "Antat.909R", "Antat.1357"), cut_line_locations = c(0, 242, 369, 693, 893, 979, 1458),
                                          cut_line_colors = c("black"), leg = "Donor VSG" , adjust_xaxis = TRUE, vert_num = 7, log_scale = FALSE, multiple = FALSE,
                                          summary = TRUE)
all_events_graph_D3_2

# summarise recomb sites per coverage, 2B

summary_guides_E9 <- all_E9 %>%
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
  filter(genotype != "NoGuide")

summary_guides_D3 <- all_D3_reseq %>%
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
  filter(genotype != "NoGuide")

E9_R1 <- read_tsv("E9_R1only_quant.tsv", col_names = c("file_info", "count"))
D3_R1 <- read_tsv("D3_R1only_quant.tsv", col_names = c("file_info", "count"))

E9_R1$file_info <- str_replace(E9_R1$file_info,"./sam_alignments/", "")
E9_R1_sum <- E9_R1 %>%  
  separate(file_info, sep = "_", into = c("mouse", "day", "genotype", "barcode", "primer", "extra" )) %>%
  select(-extra) %>%
  group_by(genotype) %>%
  summarise(total_reads = sum(count))

D3_R1$file_info <- str_replace(D3_R1$file_info,"./sam_alignments/", "")
D3_R1_sum <- D3_R1 %>%  
  separate(file_info, sep = "_", into = c("mouse", "day", "genotype", "barcode", "primer", "extra" )) %>%
  select(-extra) %>%
  group_by(genotype) %>%
  summarise(total_reads = sum(count))

E9_final <- full_join(summary_guides_E9, E9_R1_sum, by = "genotype")
D3_final <- full_join(summary_guides_D3, D3_R1_sum, by = "genotype")

E9_final_sum <- E9_final %>%
  mutate(ratio = total_reads/551398, expt = "E9") %>%
  mutate(multiplier_account = count/ratio)

D3_final_sum <- D3_final %>%
  mutate(ratio = total_reads/551398, expt = "D3") %>%
  mutate(multiplier_account = count/ratio)

cas9_ratios <- bind_rows(E9_final_sum, D3_final_sum)
cas9_ratios$genotype = str_replace(cas9_ratios$genotype, "Antat.", "")

cas9_ratios <- cas9_ratios %>%
  mutate(genotype = case_when(
    genotype == "141" ~ "guide_243",
    genotype == "300R" ~ "guide_369*",
    genotype == "592" ~ "guide_694",
    genotype == "792" ~ "guide_894",
    genotype == "909R" ~ "guide_978*",
    genotype == "1357" ~ "guide_1459",
    TRUE ~ ""
  ))

cas9_ratios <- cas9_ratios %>%
  mutate(final_log_count = log10(multiplier_account)) %>%
  qPCRr::reorder_samples("genotype", c("guide_243", "guide_369*", "guide_694", "guide_894", "guide_978*", "guide_1459"))
# 
# sum_cas9_ratios <- cas9_ratios %>%
#   group_by(genotype) %>%
#   summarise(avg = mean(final_log_count), stdev = sd(final_log_count)) #%>%
# qPCRr::reorder_samples("genotype", c("guide_243", "guide_369*", "guide_694", "guide_894", "guide_978*", "guide_1459"))
# 
# E9_norm <- cas9_ratios %>%
#   filter(expt == "E9") %>%
#   select(genotype, multiplier_account)
# D3_norm <- cas9_ratios %>%
#   filter(expt == "D3") %>%
#   select(genotype, multiplier_account)
# 
# res_aov <- aov(multiplier_account ~ genotype, data = cas9_ratios)
# summary(res_aov)
# TukeyHSD(res_aov)


ggplot(cas9_ratios) +
  geom_point(aes(x = genotype, y = final_log_count, color = expt), position = position_dodge(width = 0.5)) +
  geom_crossbar(data = sum_cas9_ratios, aes(x = genotype, y = avg, ymin = avg, ymax = avg)) +
  geom_errorbar(data = sum_cas9_ratios, aes(x = genotype, y = avg, ymin = avg - stdev, ymax = avg + stdev)) +
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
  ggplot2::labs(x = "Guides", y = "log(normalized mosaic count)", color = "Clone") +
  ggplot2::coord_cartesian(ylim = c(0, 4.5))
