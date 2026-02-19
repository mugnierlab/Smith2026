# Supplemental Figure 9A
# Mouse Infection Extended Data
# Smith 2025

# read in data
wt_consol_COs <- get_all_CO_events("./cosol_all_WT_ind_cross_over_events.csv", day_order = c("D15"),
                                   genotype_order = c("WT"), mouse_order = c("E1M7", "E1M8", "E1M10", "E2M1", "E2M2", "E2M3", "E3M6", "E3M7"),
                                   primer_order = c("1F", "2F", "3F", "4F", "0R", "1R", "2R", "3R"))

wt_consol_COs <- wt_consol_COs %>%
  filter(type_of_event %in% c("target_to_donor", "donor_to_target", "target_to_donor_edge", "donor_to_target_edge"))

mumt_consol_COs <- get_all_CO_events("./consol_all_mumt_ind_cross_over_events.csv", day_order = c("D15", "D18", "D21"),
                                     genotype_order = c("muMT"), mouse_order = c("E1M4", "E1M5", "E2M6", "E2M7", "E2M8", "E2M9", "E3M4", "M1M2",
                                                                                 "M1M3", "M1M4", "M1M5", "M1M7", "M1M8"),
                                     primer_order = c("1F", "2F", "3F", "4F", "0R", "1R", "2R", "3R"))
mumt_consol_COs <- mumt_consol_COs %>%
  filter(type_of_event %in% c("target_to_donor", "donor_to_target", "target_to_donor_edge", "donor_to_target_edge"))

all_mice <- bind_rows(wt_consol_COs, mumt_consol_COs)
all_mice <- all_mice %>%
  filter(mouse %in% c("E1M4", "E1M5", "E2M6", "E2M7", "E2M8", "E2M9", "E3M4", "E1M7", "E1M8", "E1M10", "E2M1", "E2M2", "E3M6", "E3M7"))

# normalize

wt_samples_consol_totals <- readr::read_csv("./wt_mice_consolidation_stats.txt")
mumt_samples_consol_totals <- readr::read_csv("./mumt_consolidation_stats.txt")

reads_per_sample_wt <- wt_samples_consol_totals %>%
  group_by(mouse, day, genotype) %>%
  summarise(total_reads = sum(num_consol_reads))

mumt_samples_consol_totals$genotype <- str_replace(mumt_samples_consol_totals$genotype, "muMt", "muMT")

reads_per_sample_mumt <- mumt_samples_consol_totals %>%
  group_by(mouse, day, genotype) %>%
  summarise(total_reads = sum(num_consol_reads)) %>%
  filter(day == "D15")

# norm to E1M7 consolidated
reads_per_sample_wt <- reads_per_sample_wt %>%
  mutate(normalize = 40907) %>%
  mutate(ratio = total_reads/normalize) %>%
  filter(mouse != "E2M3")

# norm to E2M7
reads_per_sample_mumt <- reads_per_sample_mumt %>%
  mutate(normalize = 2349456) %>%
  mutate(ratio = total_reads/normalize)


wt_quant <- wt_consol_COs %>%
  ungroup() %>%
  group_by(mouse, day, genotype) %>%
  summarise(total = n()) %>%
  ungroup %>%
  dplyr::add_row(mouse = "E3M6", day = "D15", genotype = "WT", total = 0)

muMT_quant <- mumt_consol_COs %>%
  ungroup() %>%
  group_by(mouse, day, genotype) %>%
  summarise(total = n()) %>%
  filter(day == "D15") %>%
  ungroup()

wt_summary_stats <- wt_quant %>%
  full_join(reads_per_sample_wt, by = c("mouse", "day", "genotype")) %>%
  mutate(norm_total = total/ratio) %>%
  filter(mouse != "E2M3")

mumt_summary_stats <- muMT_quant %>%
  full_join(reads_per_sample_mumt, by = c("mouse", "day", "genotype")) %>%
  mutate(norm_total = total/ratio)

total_sum <- bind_rows(wt_summary_stats, mumt_summary_stats)
total_summary <- total_sum %>%
  group_by(genotype, day) %>%
  summarise(avg = mean(norm_total), stdev = sd(norm_total))

total_sum <- total_sum %>%
  mutate(experiment = case_when(
    str_detect(mouse, "E1")~ "expt1",
    str_detect(mouse, "E2")~ "expt2",
    str_detect(mouse, "E3")~ "expt3"
  ))

total_sum <- total_sum %>%
  ungroup() %>%
  group_by(mouse, day, genotype, experiment) %>%
  summarise(all = sum(norm_total)) %>%
  filter(day == "D15")

total_summary <- total_sum %>%
  ungroup() %>%
  group_by(genotype, day) %>%
  summarise(avg = mean(all), stdev = sd(all)) %>%
  filter(day == "D15")

total_sum <- qPCRr::reorder_samples(total_sum, "genotype", c("WT", "muMT"))
total_summary <- qPCRr::reorder_samples(total_summary, "genotype", c("WT", "muMT"))

# graph plot
ggplot(data = total_sum) +
  geom_jitter(aes(x = day, y = all, color = experiment)) +
  geom_crossbar(data = total_summary, aes(x = day, y = avg, ymin = avg, ymax = avg)) +
  geom_errorbar(data = total_summary, aes(x = day, y = avg, ymin = avg - stdev, ymax = avg + stdev)) +
  facet_wrap(~ genotype, ncol = 2, scales = "free") +
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
