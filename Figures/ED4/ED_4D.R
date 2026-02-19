# Extended Figure 4D
# DNA double stranded breaks trigger mosaic VSG formation if homology is available
# Smith 2026

clones_inserts <- readr::read_csv("./insertion_length_all_clones.csv")

clones_inserts <- qPCRr::reorder_samples(clones_inserts, "len_desc", c("near", "far"))

clones_inserts <- clones_inserts %>%
  filter(qual == "high")

ggplot(clones_inserts) +
  geom_histogram(aes(x = insert_len), binwidth = 25) +
  facet_grid(~len_desc, scales = "free_x", space = "free_x") +
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
  ggplot2::labs(x = "Insert Length", y = "Count")

summary_stats <- clones_inserts %>%
  filter(insert_len < 300) %>%
  summarise(avg = mean(insert_len), med = median(insert_len))
