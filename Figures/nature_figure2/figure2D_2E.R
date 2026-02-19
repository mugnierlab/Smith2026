# Figure 2D & 2E
# Nucleotide changes observed in mosaic VSGs arise from templated insertions
# Smith 2026

# figure 2D
# read in and clean data, modified version of data read in to deal with alt donors available
# re-classify all VSG-2986 family members into the same purple group unless distinct at that zone

indv_cross_mechanism <-readr::read_csv("ind_cross_over_events_mech.csv",
                                       col_names = c("mouse", "day", "genotype", "primer", "direction", "read_count",
                                                     "donor_VSG", "read_start", "event_start", "event_end", "type_of_event", "UMI"),
                                       na = c("", "NA"))
indv_cross_mechanism$donor_VSG <- stringr::str_replace(indv_cross_mechanism$donor_VSG, "Tbb1125VSG-Tb427", "Tb427")
indv_cross_mechanism$donor_VSG <- stringr::str_replace(indv_cross_mechanism$donor_VSG, "Tbb1125VSG-VSG-", "VSG-")
indv_cross_mechanism$donor_VSG <- stringr::str_replace(indv_cross_mechanism$donor_VSG, "\\.1", "")

indv_cross_mechanism <- indv_cross_mechanism %>%
  group_by(mouse, day, genotype, primer, direction, read_count, donor_VSG, read_start, event_start, event_end, type_of_event, UMI) %>%
  summarise(taly = n())

weed_out_errors <- indv_cross_mechanism %>%
  dplyr::group_by(read_count, donor_VSG) %>%
  dplyr::filter(type_of_event != "donor_to_antat") %>%
  dplyr::filter(type_of_event != "antat_to_donor") %>%
  dplyr::filter(type_of_event != "antat_to_donor_edge") %>%
  dplyr::filter(type_of_event != "donor_to_antat_edge") %>%
  dplyr::summarise(count = n()) %>%
  dplyr::filter(count > 5)

trimmed_COs <- dplyr::anti_join(indv_cross_mechanism, weed_out_errors, by = "read_count")

remove_errors <- trimmed_COs %>%
  dplyr::filter(type_of_event != "random_error") %>%
  dplyr::filter(type_of_event != "random_target") %>%
  dplyr::filter(type_of_event != "random_donor")

updated_singleCOs <- remove_errors %>%
  dplyr::mutate(target_start = read_start + event_start) %>%
  dplyr::mutate(target_end = read_start + event_end) %>%
  dplyr::mutate(avg_pos = (target_start + target_end)/2)


updated_singleCOs <- updated_singleCOs %>%
  mutate(purple = case_when(
    donor_VSG == "Tb427_000362300" ~ "VSG-2986",
    donor_VSG == "Tb427_000787800" ~ "VSG-2986",
    donor_VSG == "Tb427_000143700" ~ "VSG-2986",
    TRUE ~ donor_VSG
  ))

#remove UMIs for reads that have <2 bps which match the putative donor.
updated_singleCOs <- updated_singleCOs[ ! updated_singleCOs$UMI %in% c("ACAGACATCCACGTTGTTTGGGATT", "AAATTATTATACATGTAGTATTCGG", "TATCCGATCTATCGAAGACCTAGCA",
                                                                       "ACGGTTATCTTATATACATTGCGAA", "GGCTACTCCTTGCGCCCTTTAATAA", "TTCATAACCCCTACGAAAAATGAAC"), ]



singleCOs_test <- updated_singleCOs %>%
  ungroup() %>%
  group_by(read_count, purple) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(read_count) %>%
  summarise(count_donors = n()) %>%
  filter(count_donors >1)


singleCOs <- anti_join(updated_singleCOs, singleCOs_test, by = "read_count")

singleCOs <- singleCOs %>%
  ungroup() %>%
  dplyr::select(-donor_VSG, -taly) %>%
  rename(donor_VSG = purple) %>%
  add_row(mouse = "D2-E4", day = "DOX", genotype = "NoGuide", primer = "0R", direction = "R", read_count = 6,
          read_start = 1414, event_start = 24, event_end = 25, type_of_event = "donor_to_target", 
          UMI = "TTATATGGTTCACTACGTTAGTATA", target_start = 1438, target_end = 1438.5, avg_pos = 1438.6, donor_VSG = "Ambiguous") %>%
  add_row(mouse = "D2-E4", day = "DOX", genotype = "NoGuide", primer = "0R", direction = "R", read_count = 9,
          read_start = 1428, event_start = 135, event_end = 145, type_of_event = "donor_to_target", 
          UMI = "GCAATACTTTTCGATTGGAGCCGGT", target_start = 1563, target_end = 1573, avg_pos = 1568, donor_VSG = "Ambiguous") %>%
  add_row(mouse = "D2-E4", day = "DOX", genotype = "NoGuide", primer = "0R", direction = "R", read_count = 9,
          read_start = 1428, event_start = 103, event_end = 118, type_of_event = "target_to_donor", 
          UMI = "GCAATACTTTTCGATTGGAGCCGGT", target_start = 1531, target_end = 1546, avg_pos = 1538.5, donor_VSG = "Ambiguous") %>%
  add_row(mouse = "D2-E4", day = "DOX", genotype = "NoGuide", primer = "0R", direction = "R", read_count = 10,
          read_start = 1345, event_start = 218, event_end = 228, type_of_event = "donor_to_target", 
          UMI = "CGACCCATTATCACCATATTGATCC", target_start = 1563, target_end = 1573, avg_pos = 1568, donor_VSG = "Ambiguous") %>%
  add_row(mouse = "D2-E4", day = "DOX", genotype = "NoGuide", primer = "2R", direction = "R", read_count = 18,
          read_start = 603, event_start = 92, event_end = 94, type_of_event = "donor_to_target", 
          UMI = "ACCGTTGTTGAGGATGTGACATAAT", target_start = 695, target_end = 697, avg_pos = 696, donor_VSG = "Ambiguous") %>%
  add_row(mouse = "D2-E4", day = "DOX", genotype = "NoGuide", primer = "3F", direction = "F", read_count = 25,
          read_start = 795, event_start = 72, event_end = 74, type_of_event = "target_to_donor", 
          UMI = "TACGCACACTCCAACATTTACGACC", target_start = 867, target_end = 869, avg_pos = 868, donor_VSG = "Ambiguous") %>%
  add_row(mouse = "D2-E4", day = "DOX", genotype = "Antat.592", primer = "2F", direction = "F", read_count = 168,
          read_start = 594, event_start = 139, event_end = 140, type_of_event = "target_to_donor", 
          UMI = "TGACTGCTCTGCCAAATTGGTTTTG", target_start = 733, target_end = 734, avg_pos = 733.5, donor_VSG = "Ambiguous") %>%
  add_row(mouse = "D2-E4", day = "DOX", genotype = "Antat.592", primer = "2F", direction = "F", read_count = 168,
          read_start = 594, event_start = 139, event_end = 140, type_of_event = "target_to_donor", 
          UMI = "TGACTGCTCTGCCAAATTGGTTTTG", target_start = 733, target_end = 734, avg_pos = 733.5, donor_VSG = "Ambiguous") %>%
  add_row(mouse = "D2-E8", day = "DOX", genotype = "Antat.592", primer = "1R", direction = "R", read_count = 290,
          read_start = 1186, event_start = 45, event_end = 47, type_of_event = "donor_to_target", 
          UMI = "TCTAACGCTTGGCATACTATTAACC", target_start = 1231, target_end = 1233, avg_pos = 1232, donor_VSG = "Ambiguous") %>%
  add_row(mouse = "D2-E8", day = "DOX", genotype = "Antat.592", primer = "1R", direction = "R", read_count = 290,
          read_start = 1186, event_start = 30, event_end = 30, type_of_event = "target_to_donor", 
          UMI = "TCTAACGCTTGGCATACTATTAACC", target_start = 1216, target_end = 1216, avg_pos = 1216, donor_VSG = "Ambiguous") %>%
  add_row(mouse = "D2-E8", day = "DOX", genotype = "NoGuide", primer = "0R", direction = "R", read_count = 320,
          read_start = 1421, event_start = 126, event_end = 152, type_of_event = "donor_to_target", 
          UMI = "TAAAATGAACTGCGATGATATCAAG", target_start = 1547, target_end = 1573, avg_pos = 1560, donor_VSG = "Ambiguous") %>%
  add_row(mouse = "D2-E4", day = "DOX", genotype = "Antat.592", primer = "0R", direction = "R", read_count = 334,
          read_start = 1403, event_start = 144, event_end = 170, type_of_event = "donor_to_target", 
          UMI = "GGAGTCTACAACGCTGCTGTCTGAC", target_start = 1547, target_end = 1573, avg_pos = 1560, donor_VSG = "Ambiguous")

singleCOs <- qPCRr::reorder_samples(singleCOs, "genotype", c("NoGuide", "Antat.592"))
singleCOs_E4 <- singleCOs %>%
  filter(mouse == "D2-E4")
singleCOs_E8 <- singleCOs %>%
  filter(mouse == "D2-E8")

mechan_graph <- graph_all_events(singleCOs_E4, wrap_by = "genotype",
                                 color_var = "donor_VSG", leg = "Donor VSGs", multiple = FALSE,
                                 color_scheme = TRUE, adjust_xaxis = TRUE) +
  ggplot2::geom_vline(aes(xintercept = c(86))) +
  ggplot2::geom_vline(aes(xintercept = c(1598))) +
  ggplot2::geom_vline(aes(xintercept = c(1079))) +
  ggplot2::geom_vline(aes(xintercept = c(694)))

mechan_graph_supp <- graph_all_events(singleCOs_E8, wrap_by = "genotype",
                                      color_var = "donor_VSG", leg = "Donor VSGs", multiple = FALSE,
                                      color_scheme = TRUE, adjust_xaxis = TRUE) +
  ggplot2::geom_vline(aes(xintercept = c(86))) +
  ggplot2::geom_vline(aes(xintercept = c(1598))) +
  ggplot2::geom_vline(aes(xintercept = c(1079))) +
  ggplot2::geom_vline(aes(xintercept = c(694)))


# FIgure 2E quantify 

D2_R1 <- read_tsv("D2_R1only_quant.tsv", col_names = c("file_info", "count"))

D2_R1$file_info <- str_replace(D2_R1$file_info,"./sam_alignments/", "")
D2_R1_sum <- D2_R1 %>%  
  separate(file_info, sep = "_", into = c("mouse", "day", "genotype", "barcode", "primer", "extra" )) %>%
  select(-extra) %>%
  group_by(mouse, genotype) %>%
  summarise(total_reads = sum(count))

summary_guides_D2_E4 <- singleCOs_E4 %>%
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


summary_guides_D2_E8 <- singleCOs_E8 %>%
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

quant_Norm <- as_tibble(data.frame(sample = c("D2-E4", "D2-E8", "E9", "D3"), recombinations = c(383, 805, 5263, 3983), coverage = c(817740, 832408, 707511, 903903),
                                      expt_type = c("SM", "SM", "EATRO", "EATRO"), norm = c(551398, 551398, 551398, 551398)))


quant_Norm <- qPCRr::reorder_samples(quant_Norm, "expt_type", c("EATRO", "SM"))

quant_Norm <- quant_Norm %>%
  mutate(multiplier = coverage/norm) %>%
  mutate(recomb_multiplier = recombinations/multiplier) %>%
  mutate(log_recomb = log10(recomb_multiplier))

# quant_Norm_sum <- quant_Norm %>%
#   group_by(expt_type) %>%
#   summarise(avg = mean(log_recomb), stdev = sd(log_recomb))
# 
# 
# sum <- aov(log_recomb ~expt_type, data = quant_Norm)
# summary(sum)
# TukeyHSD(sum)


ggplot(quant_Norm) +
  geom_point(aes(x = expt_type, y = log_recomb, color = sample), position = position_dodge(width = 0.5)) +
  geom_crossbar(data = quant_Norm_sum, aes(x = expt_type, y = avg, ymin = avg, ymax = avg)) +
  geom_errorbar(data = quant_Norm_sum, aes(x = expt_type, y = avg, ymin = avg - stdev, ymax = avg + stdev)) +
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
  coord_cartesian(ylim = c(0, 4))
