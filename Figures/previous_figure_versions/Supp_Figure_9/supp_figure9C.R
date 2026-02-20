# Supplemental Figure 9C
# Mouse Infection Extended Data
# Smith 2025

# D6 VS D15 staring VSGs

d6_results <- readr::read_csv("../2021_paperDraft/raw_data_Figure_invivo/results_VSGSeq/results_D6.csv",
                              col_names = c("mouse", "day", "genotype", "VSG", "percent"))

d15_results <- readr::read_csv("../2021_paperDraft/raw_data_Figure_invivo/results_VSGSeq/results_D15.csv",
                               col_names = c("mouse", "day", "genotype", "VSG", "percent"))
d15_results_mouse <- setdiff(d15_results, d6_results, by = "mouse")

results <- bind_rows(d6_results, d15_results_mouse) %>%
  filter(mouse != "E2M4") %>%
  filter(mouse != "E2M5") %>%
  filter(mouse != "E3M10") %>%
  filter(mouse != "E3M3") %>%
  filter(mouse != "E3M8") %>%
  filter(mouse != "E3M9") %>%
  filter(mouse != "EMM2") %>%
  filter(mouse != "EMM3") %>%
  filter(mouse != "EMM4") %>%
  filter(mouse != "EMM5") %>%
  filter(mouse != "EMM6") %>%
  filter(mouse != "EMM7") %>%
  filter(mouse != "EMM8") %>%
  filter(mouse != "EMM9") %>%
  filter(mouse != "EMM10") %>%
  filter(mouse != "E2M3")



other_group <- results %>%
  group_by(mouse, day, genotype) %>%
  summarise(total = sum(percent)) %>%
  mutate(VSG = "other") %>%
  mutate(percent = 100-total)

final_results <- bind_rows(results, other_group) %>%
  select(-total)
final_results <- qPCRr::reorder_samples(final_results, "day", c("D6", "D15"))
final_results <- qPCRr::reorder_samples(final_results, "genotype", c("WT", "muMT"))

ggplot(final_results) +
  geom_bar(aes(x = interaction(mouse, genotype), y = percent, color = VSG, fill = VSG ), position = "stack", stat = "identity",
           width = 0.5) +
  facet_wrap(~day)+
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
  scale_color_manual(values = c("deepskyblue", "darkgreen", "grey80"),
                     limits = c("VSG-4156", "VSG-421", "other"),
                     na.value = "grey80") +
  scale_fill_manual(values = c("deepskyblue", "darkgreen", "grey80"),
                    limits = c("VSG-4156", "VSG-421", "other"),
                    na.value = "grey80")
