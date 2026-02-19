# Figure 5C
# Mutations within VSG hypervariable region shift antibody binding
# Smith 2026

# read in data
facs_samples <- readr::read_csv("./mouse_final_calculations.csv", col_names = c("fcs", "staining", "cell", "mouse_antibody"))

controls <- facs_samples %>%
  filter(cell == "AnTat")

expts <- facs_samples %>%
  filter(cell != "AnTat")

controls_sum <- controls %>%
  group_by(mouse_antibody) %>%
  summarise(avg = mean(as.numeric(staining)), stdev = sd(as.numeric(staining)))

norm_controls <- full_join(controls, controls_sum, by = "mouse_antibody") %>%
  mutate(norm_stain = as.numeric(staining)/avg*100)

norm_controls_sum <- norm_controls %>%
  group_by(mouse_antibody, cell) %>%
  summarise(avg_final = mean(norm_stain), stdev_final = sd(norm_stain))

norm_expt <- full_join(expts, controls_sum, by = "mouse_antibody") %>%
  mutate(norm_stain = as.numeric(staining)/avg*100) 

norm_expts_sum <- norm_expt %>%
  ungroup() %>%
  group_by(mouse_antibody, cell) %>%
  summarise(avg_final = mean(norm_stain), stdev_final = sd(norm_stain))

# combined controls individual points
indiv_cntrl_points <- norm_controls %>%
  select(cell, norm_stain, mouse_antibody) %>%
  group_by(mouse_antibody, cell) %>%
  summarise(avg_final = mean(norm_stain), stdev_final = sd(norm_stain))

indiv_expt_points <- norm_expt %>%
  select(cell, norm_stain, mouse_antibody) %>%
  group_by(mouse_antibody, cell) %>%
  summarise(avg_final = mean(norm_stain), stdev_final = sd(norm_stain))

# combine controls and expts together
final_sum <- bind_rows(norm_expts_sum, norm_controls_sum)
final_points <- bind_rows(norm_controls, norm_expt)

final_sum <- final_sum %>%
  filter(cell != "VSG-221")

final_points <- final_points %>%
  filter(cell != "VSG-221")

final_points <- qPCRr::reorder_samples(final_points, "cell", c("AnTat", "VSG-228", "D6", "C10", "C3", "F7"))
final_points <- qPCRr::reorder_samples(final_points, "mouse_antibody", c("M1", "M2", "M3", "M4"))

write_csv(final_points, "facs_points.csv")

avgs_sum <- final_sum %>%
  ungroup() %>%
  group_by(cell) %>%
  summarise(total_avg = mean(avg_final), total_stdev = sd(avg_final))

avgs_sum <- qPCRr::reorder_samples(avgs_sum, "cell", c("AnTat", "VSG-228", "D6", "C10", "C3", "F7"))

mouse_color_vector = c("M1" = "darkorchid3", "M2" = "dodgerblue3", "M3" = "darkslategray4", "M4" = "palevioletred2")

#plot data <- this was used in the final figure for the paper!
ggplot(data = final_points) +
  geom_jitter(aes(x = cell, y = norm_stain, color = mouse_antibody)) +
  geom_crossbar(data = avgs_sum, aes(x = cell, y = total_avg, ymin = total_avg, ymax = total_avg)) +
  geom_errorbar(data = avgs_sum, aes(x = cell, y = total_avg, ymin = total_avg - total_stdev, ymax = total_avg + total_stdev)) +
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
  ggplot2::labs(x = "Clone", y = "Percent AnTat staining") +
  scale_y_continuous(limits = c(0, 200), breaks = seq(0,200, by=25)) +
  scale_color_manual(values = mouse_color_vector,
                     limits = c("M1", "M2", "M3", "M4"),
                     na.value = "gray")

# one way ANOVA
model <- aov(norm_stain ~ cell, data = final_points)
summary(model)
TukeyHSD(model)

