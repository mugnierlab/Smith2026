# Figure 4 A
# AnTat1.1 mosaic VSGs form in vivo and reside within tissues of WT mice
# Smith 2026

# read in data
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

# plot graph
all_mice_graph <- graph_all_events(all_mice, color_var = "donor_VSG", wrap_by = "genotype", cut_line = FALSE,
                                   leg = "Donor VSG", adjust_xaxis = TRUE, multiple = FALSE, vert_num = 7, log_scale = FALSE,
                                   mouse_plot = TRUE, rep_pic = TRUE) +
  ggplot2::geom_vline(aes(xintercept = c(86))) +
  ggplot2::geom_vline(aes(xintercept = c(1598))) +
  ggplot2::geom_vline(aes(xintercept = c(1079)))
