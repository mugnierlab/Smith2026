# Figure 5B - additional clones added
# Mutations within VSG hypervariable region shift antibody binding
# Smith 2026

# read in data
facs_samples <- readr::read_csv("./single_clones_combo.csv", col_names = c("CloneCode", "staining", "expt", "clone"))
facs_samples_round2 <- readr::read_csv("./select_rabbit_redos.csv", col_names = c("cntrl", "staining", "CloneCode", "figure_code", "donor_VSG"))

controls <- facs_samples %>%
  filter(clone == "control")

controls_round2 <- facs_samples_round2 %>%
  filter(figure_code == "Cntrl")

expts <- facs_samples %>%
  filter(clone != "control")

expts_round2 <- facs_samples_round2 %>%
  filter(figure_code != "Cntrl")


controls_sum <- controls %>%
  group_by(CloneCode) %>%
  summarise(avg = mean(staining), stdev = sd(staining))

controls_sum_round2 <- controls_round2 %>%
  group_by(cntrl) %>%
  summarise(avg = mean(staining), stdev = sd(staining))

norm_controls <- full_join(controls, controls_sum, by = "CloneCode") %>%
  mutate(norm_stain = staining/avg*100) %>%
  mutate(old_clone_name = "blagh")

norm_controls_round2 <- full_join(controls_round2, controls_sum_round2, by = "cntrl") %>%
  mutate(norm_stain = staining/avg*100)

norm_controls_sum <- norm_controls %>%
  summarise(avg_final = mean(norm_stain), stdev_final = sd(norm_stain)) %>%
  mutate(clone = "control") %>%
  mutate(old_clone_name = "blagh")

norm_controls_sum_round2 <- norm_controls_round2 %>%
  summarise(avg_final = mean(norm_stain), stdev_final = sd(norm_stain)) %>%
  mutate(clone = "control")

expts_sum <- expts %>%
  separate(CloneCode,sep = "-", c("CloneCode", "old_clone_name"))

norm_expt <- full_join(expts_sum, controls_sum, by = "CloneCode") %>%
  mutate(norm_stain = staining/avg*100) 

norm_expt_round2 <- full_join(expts_round2, controls_sum_round2, by = "cntrl") %>%
  mutate(norm_stain = staining/avg*100)

norm_expts_sum <- norm_expt %>%
  ungroup() %>%
  group_by(clone) %>%
  summarise(avg_final = mean(norm_stain), stdev_final = sd(norm_stain))

norm_expts_sum_round2 <- norm_expt_round2 %>%
  ungroup() %>%
  group_by(figure_code) %>%
  summarise(avg_final = mean(norm_stain), stdev_final = sd(norm_stain))

# combined controls individual points
indiv_cntrl_points_round1 <- norm_controls %>%
  select(clone, norm_stain, CloneCode)

indiv_cntrl_points_round2 <- norm_controls_round2 %>%
  select(cntrl, norm_stain) %>%
  rename(clone = cntrl) %>%
  mutate(CloneCode = "Expt2")

indiv_cntrl_points <- bind_rows(indiv_cntrl_points_round1, indiv_cntrl_points_round2) %>%
  mutate(figure_code = "M0") %>%
  select(-clone)

# combined controls summaries
both_cntrls_summary <- indiv_cntrl_points %>%
  group_by(figure_code) %>%
  summarise(avg_final = mean(norm_stain), stdev_final = sd(norm_stain))

# combined expts individual points

indiv_expt_points_round1 <- norm_expt %>%
  filter(old_clone_name != "H8.m") %>%
  filter(old_clone_name != "H8.228") %>%
  mutate(old_clone_name = case_when(
    old_clone_name == "B11" ~ "M3",
    old_clone_name == "A6" ~ "M2",
    old_clone_name == "A10" ~ "M1",
    old_clone_name == "C10" ~ "M4",
    TRUE ~ "None"
  )) %>%
  rename(figure_code = old_clone_name) %>%
  select(figure_code, norm_stain, CloneCode)

indiv_expt_points_round2 <- norm_expt_round2 %>%
  select(figure_code, norm_stain) %>%
  mutate(CloneCode = "Expt2")

indiv_expt_points <- bind_rows(indiv_expt_points_round1, indiv_expt_points_round2)

# combined expts summary

summary_round1 <- norm_expts_sum %>%
  filter(clone != "C2") %>%
  mutate(clone = case_when(
    clone == "C1" ~ "M1",
    clone == "C3" ~ "M2",
    clone == "C4" ~ "M3",
    clone == "C5" ~ "M4",
    TRUE ~ "None"
  )) %>%
  rename(figure_code = clone)

sum_expts_all_expts <- bind_rows(summary_round1, norm_expts_sum_round2) %>%
  filter(figure_code != "M2") %>%
  filter(figure_code != "M4")

sum_both_expts <- indiv_expt_points %>%
  filter(figure_code %in% c("M2", "M4")) %>%
  group_by(figure_code) %>%
  summarise(avg_final = mean(norm_stain), stdev_final = sd(norm_stain))

sum_expts <- bind_rows(sum_expts_all_expts, sum_both_expts)

# combine controls and expts together
final_sum <- bind_rows(sum_expts, both_cntrls_summary)
final_points <- bind_rows(indiv_cntrl_points, indiv_expt_points)

# final_points <- qPCRr::reorder_samples(final_points, "clone", c("control", "C1", "C2", "C3", "C4", "C5"))

write_csv(final_points, "facs_points.csv")

#plot data
ggplot(data = final_points) +
  geom_jitter(aes(x = figure_code, y = norm_stain)) +
  geom_crossbar(data = final_sum, aes(x = figure_code, y = avg_final, ymin = avg_final, ymax = avg_final)) +
  geom_errorbar(data = final_sum, aes(x = figure_code, y = avg_final, ymin = avg_final - stdev_final, ymax = avg_final + stdev_final)) +
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
  scale_y_continuous(limits = c(0, 130), breaks = seq(0,130, by=25))



# one way ANOVA
model <- aov(norm_stain ~ figure_code, data = final_points)
summary(model)
TukeyHSD(model)
