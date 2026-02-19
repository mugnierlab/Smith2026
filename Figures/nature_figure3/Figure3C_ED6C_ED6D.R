# Figure 3C (& Extended Data Figure 6C & 6D)
# Donor VSGs are accessible throughout the genome and selected with short, imperfect homology# Smith 2024
# Smith 2026

# read in data
all_short_trunc <- get_all_CO_events("./ind_cross_over_events_truncations.csv", day_order = c("DOX", "RAD51KO"),
                             genotype_order = c("AnTat.592.100c", "AnTat.592.200c", "AnTat.592.300c", "AnTat.592.400c", "AnTat.592.100b", "AnTat.592.200b", "AnTat.592.300b", "AnTat.592.400b", "NoGuide"), mouse_order = c("E4-C8", "Mun1"),
                             primer_order = c("1F", "2F", "3F", "4F", "0R", "1R", "2R", "3R"))

all_short_trunc_filtered <- all_short_trunc %>%
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


# for extended Data Figure 6C
graph_all_events(all_trunc_data, color_var = "donor_VSG", wrap_by = c("genotype", "mouse"), leg = "Donor VSG", adjust_xaxis = TRUE, multiple = FALSE, vert_num = 9, bin_width = 5) +
  ggplot2::geom_vline(aes(xintercept = c(86))) +
  ggplot2::geom_vline(aes(xintercept = c(1598))) +
  ggplot2::geom_vline(aes(xintercept = c(1079))) +
  ggplot2::geom_vline(aes(xintercept = c(693))) + 
  ggplot2::geom_vline(aes(xintercept = c(643))) +
  coord_cartesian(xlim = c(400, 1000))

# for Extended Data Figure 6D
subset_trunc_data <- all_trunc_data %>%
  filter(genotype %in% c("AnTat.592.Parental", "AnTat.592.100", "AnTat.592.rDNA-100-C")) %>%
  filter(donor_VSG %in% c("VSG-228")) 

subset_trunc_data$genotype <- as.character(subset_trunc_data$genotype)
subset_trunc_data <- qPCRr::reorder_samples(subset_trunc_data, "genotype", c("AnTat.592.100", "AnTat.592.rDNA-100-C"))

graph_all_events(subset_trunc_data, color_var = "donor_VSG", wrap_by = c("genotype", "mouse"), leg = "Donor VSG", multiple = TRUE, bin_width = 5) +
  ggplot2::geom_vline(aes(xintercept = c(86))) +
  ggplot2::geom_vline(aes(xintercept = c(1598))) +
  ggplot2::geom_vline(aes(xintercept = c(1079))) +
  ggplot2::geom_vline(aes(xintercept = c(693))) +
  ggplot2::geom_vline(aes(xintercept = c(677))) +
  ggplot2::geom_vline(aes(xintercept = c(777))) + 
  ggplot2::geom_vline(aes(xintercept = c(643))) +
  ggplot2::geom_vline(aes(xintercept = c(746))) + 
  coord_cartesian(xlim = c(600, 810))

# Figure 3C

all_trunc_data$mouse <- as.character(all_trunc_data$mouse)
trunc_summary <- all_trunc_data %>%
  group_by(mouse, day, genotype, donor_VSG) %>%
  filter((avg_pos >= 693 - 250) & (avg_pos <= 693 + 250)) %>%
  filter(type_of_event %in% c("target_to_donor", "donor_to_target", "target_to_donor_edge", "donor_to_target_edge")) %>%
  summarise(count = n()) %>%
  filter(genotype != "NoGuide.Parental") %>%
  filter(donor_VSG %in% c("VSG-2986", "VSG-228"))


sum_stats <- read_tsv("output_truncations_long.tsv", col_names = c("file_info", "count"))
sum_stats_second <- read_tsv("output_truncations.tsv", col_names = c("file_info", "count"))
total_sum_stats <- bind_rows(sum_stats, sum_stats_second)
total_sum_stats$file_info <- str_replace(total_sum_stats$file_info,"./sam_alignments/", "")
final_sum <- total_sum_stats %>%  
  separate(file_info, sep = "_", into = c("mouse", "day", "genotype", "barcode", "primer", "extra" )) %>%
  select(-extra) %>%
  mutate(mouse = case_when(
    barcode == "CAGAGAGG+CTAAGCCT" ~ "E4-C8b",
    barcode == "CTCTCTAC+AAGGAGTA" ~ "E4-C8b",
    barcode == "TCCTGAGC+AGAGTAGA" ~ "E4-C8a",
    barcode == "AGGCAGAA+TATCCTCT" ~ "E4-C8a",
    genotype %in% c("AnTat.592.100b", "AnTat.592.200b", "AnTat.592.300b", "AnTat.592.400b") ~ "E4-C8b",
    genotype %in% c("AnTat.592.100c", "AnTat.592.200c", "AnTat.592.300c", "AnTat.592.400c") ~ "E4-C8a",
    TRUE ~ mouse
  )) 

final_sum$genotype <- str_replace(final_sum$genotype, "b", "")
final_sum$genotype <- str_replace(final_sum$genotype, "c", "")


final_summary <- final_sum %>%
  group_by(mouse, genotype) %>%
  summarise(total_reads = sum(count)) %>%
  filter(genotype != "Antat.592") %>%
  filter(genotype != "AnTat.592.200" | mouse != "E4-C8b")


trunc_summary <- trunc_summary %>%
  filter(genotype != "AnTat.592.200" | mouse != "E4-C8b")


trunc_quant_final <- full_join(trunc_summary, final_summary, by = c("mouse", "genotype")) %>%
  mutate(ratio = count/(total_reads/420072)) %>%
  ungroup %>%
  add_row(mouse = "E4-C8a", day = "DOX", genotype = "AnTat.592.Parental", donor_VSG = "VSG-228", count = 0, total_reads = 1, ratio = 0)
trunc_quant_final <- qPCRr::reorder_samples(trunc_quant_final, "genotype", c("AnTat.592.Parental", "AnTat.592.rDNA", "AnTat.592.rDNA-noUTR-C", "AnTat.592.rDNA-500-C", "AnTat.592.400", "AnTat.592.300", "AnTat.592.200", "AnTat.592.100", "AnTat.592.rDNA-100-C"))
trunc_quant_final <- qPCRr::reorder_samples(trunc_quant_final, "donor_VSG", c("VSG-2986", "VSG-228"))

trunc_ratios_sum <- trunc_quant_final %>%
  ungroup() %>%
  group_by(genotype, donor_VSG) %>%
  summarise(avg = mean(ratio), stdev = sd(ratio))

write_csv(trunc_quant_final, "trunc_quant.csv")

ggplot(trunc_quant_final) +
  #ggplot() +
  geom_point(aes(x = genotype, y = ratio, color = mouse), position = position_dodge(width = 0.5)) +
  facet_wrap(~donor_VSG) +
  geom_crossbar(data = trunc_ratios_sum, aes(x = genotype, y = avg, ymin = avg, ymax = avg)) +
  geom_errorbar(data = trunc_ratios_sum, aes(x = genotype, y = avg, ymin = avg - stdev, ymax = avg + stdev)) +
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
  ggplot2::labs(x = "Insert Location", y = "log(normalized mosaic count)", color = "Clone") +
  ggplot2::coord_cartesian(ylim = c(0, 800))
# 
# res_aov <- aov(ratio ~ interaction(genotype, donor_VSG), data = trunc_quant_final)
# summary(res_aov)
# TukeyHSD(res_aov)



