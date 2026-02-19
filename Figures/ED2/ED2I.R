# Extended Data Figure 2I
# Cas9 system and sgRNA design and analysis
# Smith 2026
# Author: Erin Kennedy

wb_staining_puro <- tibble(c("E9", "E9", "E9", "E9", "E9", "E9", "D3", "D3", "D3", "D3", "E9", "D3", "E9", "E9", "E9"), c("g_369R","g_694", "g_894", "g_978R", "g_1459", "g_243", "g_1459", "g_369R", "g_694", "g_694", "g_1459", "g_1459", "g_369R", "g_694", "g_978R"), c("E11", "G1", "C1", "C7", "E3", "A5", "D10", "G5", "G7", "B10", "F12", "G2", "F6", "C7", "E1"),
                           c(9.43131429, 4.41325263, 2.03172136, 4.02404671, 3.32710262, 0.76643789, 3.09235599, 9.42579452, 4.73416215, 3.30417371, 3.73332776, 10.7157735, 37.0773051, 5.29851097, 9.65595945), .name_repair = ~ c("clone", "cut_site", "subclone", "normalized_intensity"))

wb_staining_puro2 <- wb_staining_puro %>%
  dplyr::ungroup() %>%
  dplyr::mutate_at("cut_site", as.factor) %>%
  dplyr::mutate(cut_site = forcats::fct_relevel(cut_site, c("g_243", "g_369R", "g_694", "g_894", "g_978R", "g_1459")))


summary_wb <- wb_staining_puro2 %>%
  group_by(cut_site) %>%
  summarise(avg = mean(normalized_intensity), stdev = sd(normalized_intensity)) %>%
  mutate(stdev = if_else(is.na(stdev), 0, stdev))

ggplot(wb_staining_puro2) +
  geom_crossbar(data = summary_wb, aes(x = cut_site, y = avg, ymin = avg, ymax = avg)) +
  geom_errorbar(data = summary_wb, aes(x = cut_site, y = avg, ymin = avg - stdev, ymax = avg + stdev)) +
  geom_jitter(aes(x = cut_site, y = normalized_intensity, color = clone)) +
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
  ggplot2::labs(x = "Guides", y = "normalized staining intensity") +
  geom_hline(aes(yintercept = 1)) +
  coord_cartesian(y = c(0, 11))
