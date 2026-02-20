# Supplemental Figure 3B+#C
# Outcomes of clonal guide selection
# Smith 2024

# read in the data supp 3B

death_counts <- readr::read_csv("./final_death_assays_counts.csv", col_names = c("Line_name", "Cas9_code", "Guide_pos",
                                                                                 "Guide_insertion_clone", "clone_count", "aimed_plate_number", "num_DOX_colonies", "num_DMSO_colonies"))

death_counts_clean <- death_counts %>%
  select(c(Guide_pos, Guide_insertion_clone, num_DOX_colonies, num_DMSO_colonies)) %>%
  mutate(ratio = num_DOX_colonies/num_DMSO_colonies*100)

death_counts_clean$Guide_pos <- as.character(death_counts_clean$Guide_pos)

death_counts_clean <- qPCRr::reorder_samples(death_counts_clean, "Guide_pos", c("243", "694", "1459"))

death_counts_summary <- death_counts_clean %>%
  group_by(Guide_pos) %>%
  summarise(avg = mean(ratio), stdev = sd(ratio))


res_aov <- aov(ratio ~ Guide_pos, data = death_counts_clean)
summary(res_aov)
TukeyHSD(res_aov)


#graph supp 3B 
ggplot(data = death_counts_clean) +
  geom_jitter(aes(x = Guide_pos, y = ratio, color = Guide_insertion_clone)) +
  geom_crossbar(data = death_counts_summary, aes(x = Guide_pos, y = avg, ymin = avg, ymax = avg)) +
  geom_errorbar(data = death_counts_summary, aes(x = Guide_pos, y = avg, ymin = avg - stdev, ymax = avg + stdev)) +
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
  ggplot2::labs(x = "Guide", y = "Mosaic Ratio", color = "clone_info")


# supp 3C Read in the data

death_outcomes <- readr::read_csv("./death_assay_outcomes.csv", col_names = TRUE)

death_outcomes_clean <- death_outcomes %>%
  mutate(Guide = case_when(
    Guide == 141 ~ 243,
    Guide == 592 ~ 694,
    Guide == 1357 ~ 1459,
    TRUE ~ 0
  )) %>%
  group_by(Guide, Outcome) %>%
  summarise(num_each = n())

death_outcomes_clean$Guide <- as.character(death_outcomes_clean$Guide)

death_outcomes_clean <- qPCRr::reorder_samples(death_outcomes_clean, "Guide", c("243", "694", "1459"))

# graph supp 3C
ggplot(death_outcomes_clean) +
  geom_bar(aes(x = Guide, y = num_each, fill = Outcome), position = "fill", stat = "identity") +
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
                 legend.text = ggplot2::element_text(size = 15))


