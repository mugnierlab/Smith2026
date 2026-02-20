#Supplemental Figure 4F
# Smith 2024

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

recomb_len <- clean_invivo_VSG8 %>%
  mutate(recomb_len = stop-start)

num_0_recombs <- recomb_len %>%
  group_by(read_name) %>%
  filter(recomb_len == 0) %>%
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

num_recombs <- cleaned_recomb %>%
  group_by(mouse, day, original_read_name) %>%
  summarise(num_events_perRead = n()) %>%
  group_by(mouse,day) %>%
  summarise(avg_recombs = median(num_events_perRead), stdev = sd(num_events_perRead))

num_mosaic_reads <- cleaned_recomb %>%
  group_by(mouse,day,original_read_name) %>%
  summarise(num_events_perRead = n())


#make the graph
graph_plot <- ggplot2::ggplot(cleaned_recomb, aes(recomb_midpoint)) +
  ggplot2::geom_histogram(ggplot2::aes(fill = mouse), binwidth = 50) +
  facet_wrap(~day, ncol = 4) +
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
  labs(x = "Sequence (bp)", y = "Count") +
  scale_x_continuous(breaks = seq(0, 1600, by = 400)) +
  coord_cartesian(xlim = c(0, 1594))
