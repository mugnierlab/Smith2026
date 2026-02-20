# Supplemental Figure 8C
#Smith 2025

rad51 <- read_xlsx("mumT_Rad51_parasitemia.xlsx", col_names = c("day_txt", "day", "mouse", "count"))
rad51 <- mutate(rad51, genotype = "muMT")
rad51 <- mutate(rad51, experiment = "E")
rad51<- mutate(rad51, count = case_when(count==0 ~ 10000,
                                                  TRUE ~ count ))
rad51 <- mutate(rad51, parasitemia_log = log(count,10))


rad51_graph <- graph_parasitemia(rad51, wrap=FALSE)

parasite_graph <- ggplot(rad51, aes(x = day, y = parasitemia_log, group = mouse)) + 
  geom_jitter(width = 0.5) +
  geom_line() +
  theme(axis.text.x = ggplot2::element_text(size = 20,
                                            angle = 90,
                                            hjust = 1,
                                            vjust = 0.35),
        axis.text.y = ggplot2::element_text(size = 20),
        axis.title.y = ggplot2::element_text(size = 25),
        axis.title.x = ggplot2::element_text(size = 25),
        plot.background = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(color = "black", size = 2),
        axis.ticks.length = ggplot2::unit(0.25, "cm"),
        legend.position = "right",
        strip.text.x = ggplot2::element_text(size = 20),
        legend.title = ggplot2::element_text(size = 15,
                                             face = "bold",
                                             hjust = 0.5),
        legend.key = ggplot2::element_rect(fill = "white"),
        legend.text = ggplot2::element_text(size = 15)) +
  labs(x = "Day", y = "Parasites/mL") +
  scale_x_continuous(breaks = seq(4, 22, by = 2)) #+
  scale_y_continuous(limits = c(4, 9), breaks = seq(4, 9, by = 1))
