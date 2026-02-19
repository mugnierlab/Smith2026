# Extended Data Figure 2E, 2F & 2G
# Cas9 system and sgRNA design and analysis
# Smith 2026

# graphing color set up
# read in data
color_VSG_vector = c("VSG-228" = "turquoise3", "VSG-2986" = "slateblue4", "VSG-3110" = "steelblue3", "VSG-7358" = "darkorchid1",
                     "VSG-4156" = "darkslategray4", "VSG-1128" = "lightpink2", "VSG-479" = "indianred4", "other VSGs" = "grey80")


# guide_infos <- readr::read_csv("./supplemental_guide_stats/guide_stats.csv",
guide_infos <- readr::read_csv("./supplemental_guide_stats/guide_stats.csv",
                               col_names = c("Guide", "donor", "pam_distance", "guide_distance",
                                             "seq_distance", "pam", "guide", "seq"))
guide_infos <- qPCRr::reorder_samples(guide_infos, "Guide", c("G_141", "G_300R", "G_592", "G_792",
                                                              "G_909R", "G_1357"))

guide_infos$donor <- str_replace(guide_infos$donor, "Tbb1125", "")

# 2E donor distance

seq_summary <- guide_infos %>%
  group_by(Guide) %>%
  summarise(avg = mean(seq_distance), stdev = sd(seq_distance))

ggplot(guide_infos) +
  geom_jitter(aes(x = Guide, y = seq_distance, color = donor), width = 0.25) +
  geom_crossbar(data = seq_summary, aes(x = Guide, y = avg, ymin = avg, ymax = avg)) +
  geom_errorbar(data = seq_summary, aes(x = Guide, y = avg, ymin = avg - stdev, ymax = avg + stdev)) +
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
                     na.value = "grey80")

# 2F guide distance

guide_summary <- guide_infos %>%
  group_by(Guide) %>%
  summarise(avg = mean(guide_distance), stdev = sd(guide_distance))

ggplot(guide_infos) +
  geom_jitter(aes(x = Guide, y = guide_distance, color = donor), width = 0.25) +
  geom_crossbar(data = guide_summary, aes(x = Guide, y = avg, ymin = avg, ymax = avg)) +
  geom_errorbar(data = guide_summary, aes(x = Guide, y = avg, ymin = avg - stdev, ymax = avg + stdev)) +
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
                     na.value = "grey80") +
  coord_cartesian(ylim = c(0,20))

# 2G PAM disrupted

pam_summary <- guide_infos %>%
  group_by(Guide) %>%
  summarise(avg = mean(pam_distance), stdev = sd(pam_distance))

pam_presence <- guide_infos %>%
  mutate(pam_intact = case_when(
    str_detect(guide_infos$pam, "[AGCT]GG") ~ "PAM",
    str_detect(guide_infos$Guide, "R") & str_detect(guide_infos$pam, "CC[ACGT]") ~ "PAM",
    TRUE ~ "not"
  )) %>%
  filter(pam_intact == "not")

ggplot(pam_presence) +
  geom_bar(aes(x = Guide, fill = donor)) +
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
  scale_fill_manual(values = color_VSG_vector,
                    limits = c("VSG-228", "VSG-2986", "VSG-3110", "VSG-7358",
                               "VSG-4156", "VSG-1128", "VSG-479", "other VSGs"),
                    na.value = "grey80")
