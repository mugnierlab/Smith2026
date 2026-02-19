# Extended Data Figure 4E & 4A
#Smith 2026

# plot ED 4E
mumtVSG8 <- readxl::read_xlsx("parasitemia_for_graphing.xlsx", col_names = TRUE)

mumtVSG8<- mutate(mumtVSG8, parasitemia = case_when(total_parasites==0 ~ 10000,
                                        TRUE ~ total_parasites ))
mumtVSG8 <- mutate(mumtVSG8, parasitemia_log = log(parasitemia,10))


parasite_graph <- ggplot(mumtVSG8, aes(x = Day, y = parasitemia_log, group = mouse)) + 
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


#plot ED 4A
wt_D3 <- readxl::read_xlsx("WT_D3_cycling_through_mice_parasitemia.xlsx", col_names = TRUE)

wt_D3<- mutate(wt_D3, parasitemia = case_when(parasitemia==0 ~ 10000,
                                                    TRUE ~ parasitemia ))
wt_D3 <- mutate(wt_D3, parasitemia_log = log(parasitemia,10))


parasite_graph <- ggplot(wt_D3, aes(x = day, y = parasitemia_log, group = mouse)) + 
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
