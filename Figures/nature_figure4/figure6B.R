# Figure 6 B
# AnTat1.1 mosaic VSGs form in vivo and reside within tissues of WT mice
# Smith 2024

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

# obtain N and C termini

#N-vs C-terminus # of donors

wt_quant_by_domain <- wt_consol_COs %>%
  ungroup() %>%
  mutate(domain = case_when(
    avg_pos <= 1080 ~ "N-term",
    avg_pos > 1080 ~ "C-term",
    TRUE ~ ""
  )) %>%
  group_by(mouse, day, genotype, domain, donor_VSG) %>%
  summarise(total = n()) %>%
  ungroup() %>%
  group_by(mouse, day, genotype, domain) %>%
  summarise(total = n()) %>%
  ungroup() %>%
  dplyr::add_row(mouse = "E3M6", day = "D15", genotype = "WT", domain = "N-term", total = 0) %>%
  dplyr::add_row(mouse = "E3M6", day = "D15", genotype = "WT", domain = "C-term", total = 0) %>%
  dplyr::add_row(mouse = "E1M8", day = "D15", genotype = "WT", domain = "N-term", total = 0) %>%
  dplyr::add_row(mouse = "E1M7", day = "D15", genotype = "WT", domain = "N-term", total = 0) %>%
  dplyr::add_row(mouse = "E1M10", day = "D15", genotype = "WT", domain = "N-term", total = 0) %>%
  dplyr::add_row(mouse = "E2M1", day = "D15", genotype = "WT", domain = "N-term", total = 0) %>%
  dplyr::add_row(mouse = "E2M2", day = "D15", genotype = "WT", domain = "N-term", total = 0)

muMT_quant_by_domain <- mumt_consol_COs %>%
  ungroup() %>%
  mutate(domain = case_when(
    avg_pos <= 1080 ~ "N-term",
    avg_pos > 1080 ~ "C-term",
    TRUE ~ ""
  )) %>%
  group_by(mouse, day, genotype, domain, donor_VSG) %>%
  summarise(total = n()) %>%
  ungroup() %>%
  group_by(mouse, day, genotype, domain) %>%
  summarise(total = n()) %>%
  ungroup() %>%
  dplyr::add_row(mouse = "E2M7", day = "D15", genotype = "muMT", domain = "N-term", total = 0)

quant_by_domain <- bind_rows(wt_quant_by_domain, muMT_quant_by_domain) %>%
  filter(mouse %in% c("E1M4", "E1M5", "E2M6", "E2M7", "E2M8", "E2M9", "E3M4", "E1M7", "E1M8", "E1M10", "E2M1", "E2M2", "E3M6", "E3M7"))


quant_by_domain <- qPCRr::reorder_samples(quant_by_domain, "genotype", c("WT", "muMT"))
quant_by_domain <- qPCRr::reorder_samples(quant_by_domain, "domain", c("N-term", "C-term"))


quant_by_domain <- quant_by_domain %>%
  mutate(experiment = case_when(
    str_detect(mouse, "E1")~ "expt1",
    str_detect(mouse, "E2")~ "expt2",
    str_detect(mouse, "E3")~ "expt3",
    str_detect(mouse, "M1M")~ "expt4",
  ))

# write_csv(muMT_quant_by_domain, "muMT_summary_stats.csv")

# write_csv(wt_quant_by_domain, "wt_summary_stats.csv")

sum_quant_domain <- quant_by_domain %>%
  ungroup() %>%
  group_by(day, genotype, domain) %>%
  summarise(avg = mean(total), stdev = sd(total))

# graph data
ggplot(data = quant_by_domain) +
  geom_jitter(aes(x = domain, y = total, color = experiment)) +
  facet_wrap(~genotype, scales = "free") +
  geom_crossbar(data = sum_quant_domain, aes(x = domain, y = avg, ymin = avg, ymax = avg)) +
  geom_errorbar(data = sum_quant_domain, aes(x = domain, y = avg, ymin = avg - stdev, ymax = avg + stdev)) +
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
  ggplot2::labs(x = "Domain", y = "Num Donor VSGs", color = "Experiment")

# stats
# nonparametric

stats::shapiro.test(quant_by_domain$total)


vsg_count_wt_N <- quant_by_domain %>%
  filter(genotype == "WT") %>%
  filter(domain == "N-term")
vsg_count_wt_N$mouse <- as_factor(vsg_count_wt_N$mouse)

vsg_count_wt_C <- quant_by_domain %>%
  filter(genotype == "WT") %>%
  filter(domain == "C-term")
vsg_count_wt_C$mouse <- as_factor(vsg_count_wt_C$mouse)

wilcox.test(vsg_count_wt_N$total, vsg_count_wt_C$total, paired = TRUE)


vsg_count_muMT_N <- quant_by_domain %>%
  filter(genotype == "muMT") %>%
  filter(domain == "N-term")
vsg_count_muMT_N$mouse <- as_factor(vsg_count_muMT_N$mouse)

vsg_count_muMT_C <- quant_by_domain %>%
  filter(genotype == "muMT") %>%
  filter(domain == "C-term")
vsg_count_muMT_C$mouse <- as_factor(vsg_count_muMT_C$mouse)

wilcox.test(vsg_count_muMT_N$total, vsg_count_muMT_C$total, paired = TRUE)
