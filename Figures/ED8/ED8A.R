# Extended Data Figure 8A
# Smith 2026

# read in data
qPCR_results <- read_csv("qPCR_for_graphing.csv")

qPCR_summary <- qPCR_results %>%
  group_by(Sample, target) %>%
  mutate(relative_expression = relative_expression * 100) %>%
  summarise(avg = mean(relative_expression), stdev = sd(relative_expression))

qPCR_summary <- qPCRr::reorder_samples(qPCR_summary, "Sample", c("AnTat J1339", "J1339 R51 -/- B1N2", "J1339 BRCA2 -/- B6N1",
                                                                 "Fx RAD51 dKO 2-1", "R51 2-1 +pJ1339 cl1", "R51 2-1 +pJ1339 cl2",
                                                                 "R51 2-1 +pJ1339 cl4", "water"))

qPCR_summary <- qPCRr::reorder_samples(qPCR_summary, "target", c("Cas9_1", "Cas9_2", "RAD51", "BRCA2"))


qPCR_for_plotting <- qPCR_results %>%
  mutate(relative_expression = relative_expression * 100)

qPCR_for_plotting <- qPCRr::reorder_samples(qPCR_for_plotting, "Sample", c("AnTat J1339", "J1339 R51 -/- B1N2", "J1339 BRCA2 -/- B6N1",
                                                                 "Fx RAD51 dKO 2-1", "R51 2-1 +pJ1339 cl1", "R51 2-1 +pJ1339 cl2",
                                                                 "R51 2-1 +pJ1339 cl4", "water"))

qPCR_for_plotting <- qPCRr::reorder_samples(qPCR_for_plotting, "target", c("Cas9_1", "Cas9_2", "RAD51", "BRCA2"))

# graph the data
ggplot(qPCR_for_plotting) +
  geom_point(aes(x = Sample, y = relative_expression)) +
  facet_wrap(~target) +
  geom_crossbar(data = qPCR_summary,
                ggplot2::aes(x = Sample, y = avg, ymin = avg, ymax = avg)) + 
  geom_errorbar(data = qPCR_summary, aes(x = Sample, y = avg,
                                         ymin = avg - stdev,
                                         ymax = avg + stdev)) +
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
