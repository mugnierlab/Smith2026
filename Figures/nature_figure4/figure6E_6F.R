# Figure 6 E & F
# AnTat1.1 mosaic VSGs form in vivo and reside within tissues of WT mice
# Smith 2024

# read in the data 
wt_consol_COs <- get_all_CO_events("./cosol_all_WT_ind_cross_over_events.csv", day_order = c("D15"),
                                   genotype_order = c("WT"), mouse_order = c("E1M7", "E1M8", "E1M10", "E2M1", "E2M2", "E2M3", "E3M6", "E3M7"),
                                   primer_order = c("1F", "2F", "3F", "4F", "0R", "1R", "2R", "3R"))

wt_consol_COs <- wt_consol_COs %>%
  filter(type_of_event %in% c("target_to_donor", "donor_to_target", "target_to_donor_edge", "donor_to_target_edge"))

mumt_consol_COs <- get_all_CO_events("./consol_all_mumt_ind_cross_over_events.csv", day_order = c("D15", "D18", "D21"),
                                     genotype_order = c("muMT"), mouse_order = c("E1M4", "E1M5", "E2M6", "E2M7", "E2M8", "E2M9", "E3M4", "M1M2",
                                                                                 "M1M3", "M1M4", "M1M5", "M1M7", "M1M8"),
                                     primer_order = c("1F", "2F", "3F", "4F", "0R", "1R", "2R", "3R"))
mumt_consol_COs <- mumt_consol_COs %>%
  filter(type_of_event %in% c("target_to_donor", "donor_to_target", "target_to_donor_edge", "donor_to_target_edge"))

all_mice <- bind_rows(wt_consol_COs, mumt_consol_COs)
all_mice <- all_mice %>%
  filter(mouse %in% c("E1M4", "E1M5", "E2M6", "E2M7", "E2M8", "E2M9", "E3M4", "E1M7", "E1M8", "E1M10", "E2M1", "E2M2", "E3M6", "E3M7"))


# 6E recombination site length

wt_recomb_sites <- wt_consol_COs %>%
  mutate(recomb_len = target_end - target_start)

recomb_mumt <- mumt_consol_COs %>%
  mutate(recomb_len = target_end - target_start)

recomb_mice <- bind_rows(recomb_mumt, wt_recomb_sites)

all_recomb_sites_mice <- ggplot(recomb_mice) +
  geom_histogram(aes(x = recomb_len), binwidth = 5) +
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
  ggplot2::labs(x = "Recomb Site Length", y = "Count") #+
coord_cartesian(xlim = c(0, 200))

# 6F insertion lengths

# insertion lengths ie double cross overs
num_mice_reads <- all_mice %>%
  group_by(genotype, read_count) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(genotype) %>%
  summarise(count = n())

inserts_tally_mice <- all_mice %>%
  group_by(genotype, read_count) %>%
  summarise(count = n()) %>%
  filter(count > 1) %>%
  ungroup() %>%
  group_by(genotype) %>%
  summarise(count = n())


"%notin%" <- Negate("%in%")
# E9 % donors that are within N-term which are not a family member
N_term_mice <- all_mice %>%
  filter(avg_pos <= 1079)
#WT = 2
#mumt = 5991
N_term_mice_alt <- N_term_mice %>%
  group_by(genotype, donor_VSG) %>%
  filter(donor_VSG %notin% c("VSG-4156", "VSG-7358", "VSG-3110", "VSG-228", "VSG-2986", "family"))
#WT = 0%
#muMT = 5, 0.08% 99.9%

#insert_lengths
nonInsert_nonmulti_reads_mice <- all_mice %>%
  group_by(genotype, read_count) %>%
  summarise(count = n()) %>%
  filter(count != 2)

multi_inserts_mice <- all_mice %>%
  group_by(genotype, read_count) %>%
  summarise(count = n()) %>%
  filter(count <= 2)

multi_inserts_mice <- anti_join(all_mice, multi_inserts_mice, by = c("read_count", "genotype"))

multi_inserts_even_mice <- multi_inserts_mice %>%
  group_by(genotype, read_count) %>%
  summarise(count = n()) %>%
  mutate(remainder = count %% 2) %>%
  filter(remainder == 0)

multi_inserts_odd_mice <- anti_join(multi_inserts_mice, multi_inserts_even_mice, by = c("read_count", "genotype")) %>%
  group_by(genotype, read_count) %>%
  mutate(rank = order(avg_pos)) %>%
  mutate(new_rank = case_when(
    rank == as.integer(1) & type_of_event %in% c("target_to_donor", "target_to_donor_edge") ~ rank,
    rank == as.integer(1) & type_of_event %in% c("donor_to_target", "donor_to_target_edge") ~ as.integer(0),
    TRUE ~ rank
  )) %>%
  filter(new_rank != as.integer(0))

start_with_target_mice <- multi_inserts_odd_mice %>%
  summarise(count = n()) %>%
  mutate(remainder = count %% 2) %>%
  filter(remainder > 0)

start_with_donor_mice <- anti_join(multi_inserts_odd_mice, start_with_target_mice, by = c("read_count", "genotype")) %>%
  select(-rank, -new_rank)

multi_inserts_odd_targetStart_mice <- anti_join(multi_inserts_mice, multi_inserts_even_mice, by = c("read_count", "genotype")) %>%
  anti_join(start_with_donor_mice, by = c("read_count", "genotype")) %>%
  group_by(genotype, read_count) %>%
  mutate(rank = order(avg_pos, decreasing = TRUE)) %>%
  mutate(new_rank = case_when(
    rank == as.integer(1) & type_of_event %in% c("target_to_donor", "target_to_donor_edge") ~ as.integer(0),
    TRUE ~ rank
  )) %>%
  filter(new_rank != as.integer(0)) %>%
  select(-rank, -new_rank)

# read_count 4294 has 2 cross overs.... will just remove and add back for now
insert_reads_mice <- anti_join(all_mice, nonInsert_nonmulti_reads_mice, by = c("read_count", "genotype")) %>%
  bind_rows(start_with_donor_mice, multi_inserts_odd_targetStart_mice) %>%
  filter(read_count != 14865) %>%
  filter(read_count != 16498) %>%
  filter(read_count != 35954) %>%
  filter(read_count != 35960) %>%
  filter(read_count != 16548) %>%
  filter(read_count != 16648) %>%
  filter(read_count != 34042) %>%
  filter(read_count != 407) %>%
  filter(read_count != 385) %>%
  filter(read_count != 1441) %>%
  filter(read_count != 1770) %>%
  filter(read_count != 5136) %>%
  filter(read_count != 7777) %>%
  filter(read_count != 26578) %>%
  filter(read_count != 26631) %>%
  filter(read_count != 31951)

# need to add these back!
starts_mice <- insert_reads_mice %>%
  filter(type_of_event %in% c("target_to_donor", "target_to_donor_edge"))
ends_mice <- insert_reads_mice %>%
  filter(type_of_event %in% c("donor_to_target", "donor_to_target_edge"))

starts_mice <- starts_mice %>%
  select(genotype, read_count, target_end)
ends_mice <- ends_mice %>%
  select(genotype, read_count, target_start)

distances_mice <- left_join(starts_mice, ends_mice, by = c("genotype", "read_count")) %>%
  add_row(genotype = "muMT", read_count = 14865, target_end = 1348, target_start = 1364) %>%
  add_row(genotype = "muMT", read_count = 16498, target_end = 927, target_start = 944) %>%
  add_row(genotype = "muMT", read_count = 35954, target_end = 735, target_start = 762) %>%
  add_row(genotype = "muMT", read_count = 35960, target_end = 734, target_start = 762) %>%
  add_row(genotype = "muMT", read_count = 16548, target_end = 841, target_start = 843) %>%
  add_row(genotype = "muMT", read_count = 16648, target_end = 841, target_start = 843) %>%
  add_row(genotype = "muMT", read_count = 34042, target_end = 475, target_start = 482) %>%
  add_row(genotype = "muMT", read_count = 407, target_end = 1298, target_start = 1313) %>%
  add_row(genotype = "WT", read_count = 385, target_end = 1298, target_start = 1330) %>%
  add_row(genotype = "muMT", read_count = 385, target_end = 1428, target_start = 1457) %>%
  add_row(genotype = "muMT", read_count = 1441, target_end = 1529, target_start = 1533) %>%
  add_row(genotype = "muMT", read_count = 1441, target_end = 1545, target_start = 1547) %>%
  add_row(genotype = "muMT", read_count = 1770, target_end = 1529, target_start = 1533) %>%
  add_row(genotype = "muMT", read_count = 1770, target_end = 1545, target_start = 1547) %>%
  add_row(genotype = "muMT", read_count = 5136, target_end = 1529, target_start = 1533) %>%
  add_row(genotype = "muMT", read_count = 5136, target_end = 1540, target_start = 1547) %>%
  add_row(genotype = "muMT", read_count = 7777, target_end = 1529, target_start = 1533) %>%
  add_row(genotype = "muMT", read_count = 7777, target_end = 1545, target_start = 1547) %>%
  add_row(genotype = "muMT", read_count = 26578, target_end = 1529, target_start = 1533) %>%
  add_row(genotype = "muMT", read_count = 26578, target_end = 1545, target_start = 1547) %>%
  add_row(genotype = "muMT", read_count = 26631, target_end = 1529, target_start = 1533) %>%
  add_row(genotype = "muMT", read_count = 26631, target_end = 1545, target_start = 1547) %>%
  add_row(genotype = "muMT", read_count = 31951, target_end = 1529, target_start = 1533) %>%
  add_row(genotype = "muMT", read_count = 31951, target_end = 1545, target_start = 1547) %>%
  mutate(insert_len = target_start - target_end) %>%
  # remove reads where it appears antat is inserted into the donor
  filter(insert_len > 0)

distances_sum_mice <- distances_mice %>%
  summarise(avg = mean(insert_len), stdev = sd(insert_len), med = median(insert_len))

insert_dis_mice <- ggplot(distances_mice) +
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
