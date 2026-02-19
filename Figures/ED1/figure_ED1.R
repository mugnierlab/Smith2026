# Extended Data Figure 1A & 1B
# Mosaic Recombination site & VSG-AMP-Seq Validation
# Smith 2026

# 1A read in data
control_results <- readr::read_csv("./percent_VSGs.csv",
                                   col_names = c("mouse", "day", "genotype", "VSG", "percent"))

control_results <- qPCRr::reorder_samples(control_results, "VSG", c("other", "VSG-228", "AnTat"))

# graph results
ggplot(control_results) +
  geom_bar(aes(x = interaction(mouse, genotype), y = percent, color = VSG, fill = VSG ), position = "stack", stat = "identity") +
  facet_wrap(~day)+
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
  scale_color_manual(values = c("deepskyblue", "darkgreen", "grey80"),
                     limits = c("AnTat", "VSG-228", "other"),
                     na.value = "grey80") +
  scale_fill_manual(values = c("deepskyblue", "darkgreen", "grey80"),
                    limits = c("AnTat", "VSG-228", "other"),
                    na.value = "grey80")
# 1B
all_SM_spike <- get_all_CO_events("./ind_cross_over_events_small_expts.csv", day_order = c("test"),
                                  genotype_order = c("spike"), mouse_order = c("SM-Antat-228", "SM-Antat"),
                                  primer_order = c("1F", "2F", "3F", "4F", "0R", "1R", "2R", "3R"))
all_SM_spike_noErrors <- all_SM_spike %>%
  filter(type_of_event %in% c("donor_to_target", "target_to_donor", "target_to_donor_edge", "donor_to_target_edge"))

graph_all_events(all_SM_spike_noErrors, color_var = "donor_VSG", wrap_by = "mouse", leg = "Donor VSG", adjust_xaxis = TRUE, multiple = TRUE) +
  ggplot2::geom_vline(aes(xintercept = c(86))) +
  ggplot2::geom_vline(aes(xintercept = c(1598))) +
  ggplot2::geom_vline(aes(xintercept = c(1079)))
