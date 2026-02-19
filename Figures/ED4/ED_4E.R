# Extended Data Figure 4E
# DNA breaks result in switchers if a homologous donor is not present
# Smith 2026

# Clone 1 is D3, Clone 2 is E9

#read in data:
all_E9 <- get_all_CO_events("./ind_cross_over_events_E9.csv", day_order = c("DOX", "DMSO"),
                            genotype_order = c("NoGuide", "Antat.141", "Antat.300R", "Antat.592",
                                               "Antat.792", "Antat.909R", "Antat.1357"), mouse_order = c("E9"),
                            primer_order = c("1F", "2F", "3F", "4F", "0R", "1R", "2R", "3R"))

all_D3_reseq <- get_all_CO_events("./ind_cross_over_events_D32.csv", day_order = c("DOX", "DMSO"),
                                  genotype_order = c("NoGuide", "Antat.141", "Antat.300R", "Antat.592",
                                                     "Antat.792", "Antat.909R", "Antat.1357"), mouse_order = c("D3"),
                                  primer_order = c("1F", "2F", "3F", "4F", "0R", "1R", "2R", "3R"))

all_E9 <- all_E9 %>%
  filter(type_of_event %in% c("target_to_donor", "donor_to_target", "target_to_donor_edge", "donor_to_target_edge")) %>%
  filter(day != "DMSO")

all_D3_reseq <- all_D3_reseq %>%
  filter(type_of_event %in% c("target_to_donor", "donor_to_target", "target_to_donor_edge", "donor_to_target_edge")) %>%
  filter(day != "DMSO")

# find insert lengths
# inserts for D3  

nonInsert_nonmulti_reads_D3 <- all_D3_reseq %>%
  group_by(read_count) %>%
  summarise(count = n()) %>%
  filter(count != 2)

multi_inserts_D3 <- all_D3_reseq %>%
  group_by(read_count) %>%
  summarise(count = n()) %>%
  filter(count <= 2)

multi_inserts_D3 <- anti_join(all_D3_reseq, multi_inserts_D3, by = "read_count")

multi_inserts_even_D3 <- multi_inserts_D3 %>%
  group_by(read_count) %>%
  summarise(count = n()) %>%
  mutate(remainder = count %% 2) %>%
  filter(remainder == 0)

multi_inserts_odd_D3 <- anti_join(multi_inserts_D3, multi_inserts_even_D3, by = "read_count") %>%
  group_by(read_count) %>%
  mutate(rank = order(avg_pos)) %>%
  mutate(new_rank = case_when(
    rank == as.integer(1) & type_of_event %in% c("target_to_donor", "target_to_donor_edge") ~ rank,
    rank == as.integer(1) & type_of_event %in% c("donor_to_target", "donor_to_target_edge") ~ as.integer(0),
    TRUE ~ rank
  )) %>%
  filter(new_rank != as.integer(0))
start_with_target_D3 <- multi_inserts_odd_D3 %>%
  summarise(count = n()) %>%
  mutate(remainder = count %% 2) %>%
  filter(remainder > 0)

start_with_donor_D3 <- anti_join(multi_inserts_odd_D3, start_with_target_D3, by = "read_count") %>%
  select(-rank, -new_rank)

multi_inserts_odd_targetStart_D3 <- anti_join(multi_inserts_D3, multi_inserts_even_D3, by = "read_count") %>%
  anti_join(start_with_donor_D3, by = "read_count") %>%
  group_by(read_count) %>%
  mutate(rank = order(avg_pos, decreasing = TRUE)) %>%
  mutate(new_rank = case_when(
    rank == as.integer(1) & type_of_event %in% c("target_to_donor", "target_to_donor_edge") ~ as.integer(0),
    TRUE ~ rank
  )) %>%
  filter(new_rank != as.integer(0)) %>%
  select(-rank, -new_rank)

insert_reads_D3 <- anti_join(all_D3_reseq, nonInsert_nonmulti_reads_D3, by = "read_count") %>%
  bind_rows(start_with_donor_D3, multi_inserts_odd_targetStart_D3)

starts_D3 <- insert_reads_D3 %>%
  filter(type_of_event %in% c("target_to_donor", "target_to_donor_edge"))
ends_D3 <- insert_reads_D3 %>%
  filter(type_of_event %in% c("donor_to_target", "donor_to_target_edge"))

starts_D3 <- starts_D3 %>%
  select(read_count, target_end)
ends_D3 <- ends_D3 %>%
  select(read_count, target_start)
distances_D3 <- left_join(starts_D3, ends_D3, by = c("read_count")) %>%
  mutate(insert_len = target_start - target_end) %>%
  filter(insert_len > 0)

distances_sum_D3 <- distances_D3 %>%
  summarise(avg = mean(insert_len), stdev = sd(insert_len), med = median(insert_len))

# insert_lengths E9

nonInsert_nonmulti_reads_E9 <- all_E9 %>%
  group_by(read_count) %>%
  summarise(count = n()) %>%
  filter(count != 2)

multi_inserts <- all_E9 %>%
  group_by(read_count) %>%
  summarise(count = n()) %>%
  filter(count <= 2)

multi_inserts <- anti_join(all_E9, multi_inserts, by = "read_count")

multi_inserts_even <- multi_inserts %>%
  group_by(read_count) %>%
  summarise(count = n()) %>%
  mutate(remainder = count %% 2) %>%
  filter(remainder == 0)

multi_inserts_odd <- anti_join(multi_inserts, multi_inserts_even, by = "read_count") %>%
  group_by(read_count) %>%
  mutate(rank = order(avg_pos)) %>%
  mutate(new_rank = case_when(
    rank == as.integer(1) & type_of_event %in% c("target_to_donor", "target_to_donor_edge") ~ rank,
    rank == as.integer(1) & type_of_event %in% c("donor_to_target", "donor_to_target_edge") ~ as.integer(0),
    TRUE ~ rank
  )) %>%
  filter(new_rank != as.integer(0))
start_with_target <- multi_inserts_odd %>%
  summarise(count = n()) %>%
  mutate(remainder = count %% 2) %>%
  filter(remainder > 0)
start_with_donor <- anti_join(multi_inserts_odd, start_with_target, by = "read_count") %>%
  select(-rank, -new_rank)

multi_inserts_odd_targetStart <- anti_join(multi_inserts, multi_inserts_even, by = "read_count") %>%
  anti_join(start_with_donor, by = "read_count") %>%
  group_by(read_count) %>%
  mutate(rank = order(avg_pos, decreasing = TRUE)) %>%
  mutate(new_rank = case_when(
    rank == as.integer(1) & type_of_event %in% c("target_to_donor", "target_to_donor_edge") ~ as.integer(0),
    TRUE ~ rank
  )) %>%
  filter(new_rank != as.integer(0)) %>%
  select(-rank, -new_rank)

# read_count 4294 has 2 cross overs.... will just remove and add back for now
insert_reads_E9 <- anti_join(all_E9, nonInsert_nonmulti_reads_E9, by = "read_count") %>%
  bind_rows(start_with_donor, multi_inserts_odd_targetStart) %>%
  filter(read_count != 4294)
starts_E9 <- insert_reads_E9 %>%
  filter(type_of_event %in% c("target_to_donor", "target_to_donor_edge"))
ends_E9 <- insert_reads_E9 %>%
  filter(type_of_event %in% c("donor_to_target", "donor_to_target_edge"))

starts_E9 <- starts_E9 %>%
  select(read_count, target_end)
ends_E9 <- ends_E9 %>%
  select(read_count, target_start)
distances_E9 <- left_join(starts_E9, ends_E9, by = c("read_count")) %>%
  add_row(read_count = 4294, target_end = 786, target_start = 788) %>%
  add_row(read_count = 4294, target_end = 793, target_start = 795) %>%
  mutate(insert_len = target_start - target_end) %>%
  filter(insert_len > 0)

distances_sum_E9 <- distances_E9 %>%
  summarise(avg = mean(insert_len), stdev = sd(insert_len), med = median(insert_len))

# graph inserts

insert_dis_D3 <- ggplot(distances_D3) +
  geom_histogram(aes(x = insert_len), binwidth = 25) +
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
  ggplot2::labs(x = "Insert Length", y = "Count") +
  coord_cartesian(xlim = c(0, 200))

insert_dis_E9 <- ggplot(distances_E9) +
  geom_histogram(aes(x = insert_len), binwidth = 25) +
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
  ggplot2::labs(x = "Insert Length", y = "Count") +
  coord_cartesian(xlim = c(0, 200))
