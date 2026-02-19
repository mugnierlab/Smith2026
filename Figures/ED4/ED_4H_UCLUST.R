# Extended Data Figure 4H (right, UCLUST)
# DNA breaks result in switchers if homologous donor is not present
# Smith 2026

VSG_group_counts <- tibble(c("EATRO", "EATRO", "EATRO", "Lister", "Lister", "Lister"), c("Family", "Unique","Duplicate", "Family", "Unique", "Duplicate"),
                           c(5268-996, 996, 0, 5667-902-4, 902-4, 4), .name_repair = ~ c("cell_line", "type", "count"))
VSG_group_counts <- VSG_group_counts %>%
  qPCRr::reorder_samples("type", c("Unique", "Duplicate", "Family"))


ggplot(VSG_group_counts, aes(fill = type, y = count, x = cell_line)) + 
  geom_bar(position="fill", stat="identity") +
  ggplot2::labs(x = "Cell Lines", y = "Percent", fill = "VSG Type") +
  theme(axis.text.x = ggplot2::element_text(size = 20,
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
  scale_fill_manual(values = c("grey80", "grey30", "grey60"))
