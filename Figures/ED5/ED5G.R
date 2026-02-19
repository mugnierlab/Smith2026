#Supplemental Figure 5G
# VSG-8 infections in muMT mice
# Smith 2026

# read in data
corrected_insertions_VSG8_invivo <- read_tsv("all_mos_cross_overs.tsv", col_names = c("read_name", "start", "stop", "type_of_event"))

read_name_sample <- read_tsv("read_name_sampleinfo.tsv", col_names = c("original_read_name", "sample_info", "read_name"))

full_vsg8_invivo <- full_join(read_name_sample, corrected_insertions_VSG8_invivo, by = "read_name")

full_vsg8_invivo_clean <- full_vsg8_invivo %>%
  drop_na()

clean_invivo_VSG8 <- full_vsg8_invivo_clean %>%
  filter(type_of_event %in% c("target_to_donor", "donor_to_target", "target_to_donor_edge", "donor_to_target_edge"))

error_count <- full_vsg8_invivo_clean %>%
  filter(type_of_event %in% c("random_error", "random_target", "random_donor"))

summary_num_insertions <- clean_invivo_VSG8 %>%
  group_by(read_name) %>%
  summarise(count = n())

homology_len <- clean_invivo_VSG8 %>%
  mutate(homology_len = stop-start)

num_0_recombs <- homology_len %>%
  group_by(read_name) %>%
  filter(homology_len == 0) %>%
  summarise(count = n()) %>%
  filter(count >= 3)

remove_extra_recombs <- anti_join(clean_invivo_VSG8, num_0_recombs, by = "read_name")
summary_num_insertions_clean <- remove_extra_recombs %>%
  group_by(read_name) %>%
  summarise(count = n())

summary_num_errors <- error_count %>%
  group_by(read_name) %>%
  summarise(count = n()) %>%
  filter(count >= 100)

remove_extra_recombs_errorsRemoved <- anti_join(remove_extra_recombs, summary_num_errors, by = "read_name")

cleaned_recomb <- remove_extra_recombs_errorsRemoved %>%
  mutate(recomb_midpoint = (stop + start)/2) %>%
  separate_wider_delim(sample_info, delim = "_", names = c("mouse", "day")) %>%
  reorder_samples("day", c("D6", "D13", "D15", "D18"))

# normalize data
num_mosaic_reads <- cleaned_recomb %>%
  group_by(mouse,day,original_read_name) %>%
  summarise(num_events_perRead = n()) %>%
  group_by(mouse,day) %>%
  summarise(mosaic_read_count = n())

norm_counts <- read_tsv("loading_VSG8_Sequences.tsv", col_names = c("sample_info", "count"))

norm_counts <- norm_counts %>%
  separate_wider_delim(sample_info, delim = "_", names = c("mouse", "day"))

recomb_normed <- full_join(num_mosaic_reads, norm_counts, by = c("mouse", "day"))
# norm to M4D13 

recomb_normed <- recomb_normed %>%
  mutate(norm_recombs = mosaic_read_count/(count/4247021)) %>%
  reorder_samples("day", c("D6", "D13", "D15", "D18"))

recomb_normed_sum <- recomb_normed %>%
  group_by(day) %>%
  summarise(avg = mean(norm_recombs), stdev = sd(norm_recombs))

# graph the data

ggplot(recomb_normed) +
  #ggplot() +
  geom_point(aes(x = day, y = norm_recombs, color = mouse), position = position_dodge(width = 0.5)) +
  geom_crossbar(data = recomb_normed_sum, aes(x = day, y = avg, ymin = avg, ymax = avg)) +
  geom_errorbar(data = recomb_normed_sum, aes(x = day, y = avg, ymin = avg - stdev, ymax = avg + stdev)) +
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
  ggplot2::labs(x = "Insert Location", y = "", color = "Clone") #+
# ggplot2::coord_cartesian(ylim = c(0, 4.5))

# stats
# nonparametric

stats::shapiro.test(recomb_normed$norm_recombs)


recomb_normed_D6 <- recomb_normed %>%
  filter(day == "D6")
recomb_normed_D6$day <- as_factor(recomb_normed_D6$day)

recomb_normed_late <- recomb_normed %>%
  filter(day != "D6") %>%
  filter(day != "D18")
recomb_normed_late$day <- as_factor(recomb_normed_late$day)


wilcox.test(recomb_normed_D6$norm_recombs, recomb_normed_late$norm_recombs, paired = TRUE)

