# Figure 4 C & D
# AnTat1.1 mosaic VSGs form in vivo and reside within tissues of WT mice
# Smith 2026

# read in data

full_tissues_COs_final <- get_all_CO_events("./single_COs_withCterm_cleaned.csv", day_order = c("D6", "D10", "D14"),
                                            genotype_order = c("Blood", "Brain", "Ear", "GonFat", "SubCu", "Heart", "Lung"),
                                            mouse_order = c("M30", "M36", "M4", "M6", "M2", "C3M4", "M34", "M3", "M12"),
                                            primer_order = c("1R"))




full_tissues_COs_final <- full_tissues_COs_final %>%
  filter(mouse != "C3M4") %>%
  filter(read_count != 6)

# 4C
# graph the tissues 
tissue_VSG_seq_graphs_final <- graph_all_events(full_tissues_COs_final, wrap_by = "tissues",
                                                color_var = "donor_VSG", leg = "Donor VSGs", multiple = FALSE,
                                                color_scheme = TRUE, adjust_xaxis = TRUE, by_mouse_only = TRUE) +
  ggplot2::geom_vline(aes(xintercept = c(86))) +
  ggplot2::geom_vline(aes(xintercept = c(1598))) +
  ggplot2::geom_vline(aes(xintercept = c(1079)))

# 4D
tissues_quant_by_point <- full_tissues_COs_final %>%
  mutate(tissues = case_when(
    genotype == "Blood" ~ "Blood",
    TRUE ~ "Tissues")) %>%
  group_by(day, mouse, tissues) %>%
  summarise(total = n()) %>%
  ungroup() %>%
  dplyr::add_row(day = "D6", mouse = "M12", tissues = "Blood", total = 0) %>%
  dplyr::add_row(day = "D6", mouse = "M27", tissues = "Blood", total = 0) %>%
  dplyr::add_row(day = "D6", mouse = "M32", tissues = "Blood", total = 0) %>%
  dplyr::add_row(day = "D6", mouse = "M33", tissues = "Blood", total = 0) %>%
  dplyr::add_row(day = "D6", mouse = "M4", tissues = "Blood", total = 0) %>%
  dplyr::add_row(day = "D6", mouse = "M30", tissues = "Blood", total = 0) %>%
  dplyr::add_row(day = "D6", mouse = "M36", tissues = "Blood", total = 0) %>%
  dplyr::add_row(day = "D6", mouse = "M6", tissues = "Blood", total = 0) %>%
  dplyr::add_row(day = "D6", mouse = "M3", tissues = "Blood", total = 0) %>%
  dplyr::add_row(day = "D6", mouse = "M2", tissues = "Blood", total = 0) %>%
  dplyr::add_row(day = "D6", mouse = "M34", tissues = "Blood", total = 0) %>%
  dplyr::add_row(day = "D6", mouse = "M35", tissues = "Blood", total = 0) %>%
  dplyr::add_row(day = "D6", mouse = "M27", tissues = "Tissues", total = 0) %>%
  dplyr::add_row(day = "D6", mouse = "M32", tissues = "Tissues", total = 0) %>%
  dplyr::add_row(day = "D6", mouse = "M33", tissues = "Tissues", total = 0) %>%
  dplyr::add_row(day = "D10", mouse = "M4", tissues = "Blood", total = 0) %>%
  dplyr::add_row(day = "D10", mouse = "M30", tissues = "Blood", total = 0) %>%
  dplyr::add_row(day = "D10", mouse = "M6", tissues = "Blood", total = 0) %>%
  dplyr::add_row(day = "D10", mouse = "M3", tissues = "Blood", total = 0) %>%
  dplyr::add_row(day = "D10", mouse = "M2", tissues = "Blood", total = 0) %>%
  dplyr::add_row(day = "D10", mouse = "M35", tissues = "Blood", total = 0) %>%
  dplyr::add_row(day = "D14", mouse = "M2", tissues = "Blood", total = 0) %>%
  dplyr::add_row(day = "D14", mouse = "M35", tissues = "Blood", total = 0) %>%
  dplyr::add_row(day = "D14", mouse = "M34", tissues = "Blood", total = 0) %>%
  dplyr::add_row(day = "D14", mouse = "M35", tissues = "Tissues", total = 0) %>%
  qPCRr::reorder_samples("day", c("D6", "D10", "D14"))

tissues_point_sum <- tissues_quant_by_point %>%
  group_by(day, tissues) %>%
  summarise(avg = mean(total), stdev = sd(total))

ggplot(data = tissues_quant_by_point) +
  geom_jitter(aes(x = day, y = total, color = tissues), position = position_jitterdodge(dodge.width = 1)) +
  geom_crossbar(data = tissues_point_sum, aes(x = day, y = avg, ymin = avg, ymax = avg, color = tissues), position = position_dodge(width = 1)) +
  geom_errorbar(data = tissues_point_sum, aes(x = day, y = avg, ymin = avg - stdev, ymax = avg + stdev, color = tissues), position = position_dodge(width = 1)) +
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
  ggplot2::labs(x = "Day", y = "# Recombinations", color = "Source") 

stats::shapiro.test(tissues_quant_by_point$total)

tissues_quant_by_point_D6 <- tissues_quant_by_point %>%
  filter(day == "D6")
tissues_quant_by_point_D6_blood <- tissues_quant_by_point_D6 %>%
  filter(tissues == "Blood")
tissues_quant_by_point_D6_tissues <- tissues_quant_by_point_D6 %>%
  filter(tissues == "Tissues")
tissues_quant_by_point_D6_blood$mouse <- as_factor(tissues_quant_by_point_D6_blood$mouse)
tissues_quant_by_point_D6_tissues$mouse <- as_factor(tissues_quant_by_point_D6_tissues$mouse)


wilcox.test(tissues_quant_by_point_D6_blood$total, tissues_quant_by_point_D6_tissues$total, p.adjust.methods = "BH")

tissues_quant_by_point_D10 <- tissues_quant_by_point %>%
  filter(day == "D10")
tissues_quant_by_point_D10_blood <- tissues_quant_by_point_D10 %>%
  filter(tissues == "Blood")
tissues_quant_by_point_D10_tissues <- tissues_quant_by_point_D10 %>%
  filter(tissues == "Tissues")
tissues_quant_by_point_D10_blood$mouse <- as_factor(tissues_quant_by_point_D10_blood$mouse)
tissues_quant_by_point_D10_tissues$mouse <- as_factor(tissues_quant_by_point_D10_tissues$mouse)
wilcox.test(tissues_quant_by_point_D10_blood$total, tissues_quant_by_point_D10_tissues$total, p.adjust.methods = "BH")


tissues_quant_by_point_D14 <- tissues_quant_by_point %>%
  filter(day == "D14")
tissues_quant_by_point_D14_blood <- tissues_quant_by_point_D14 %>%
  filter(tissues == "Blood")
tissues_quant_by_point_D14_tissues <- tissues_quant_by_point_D14 %>%
  filter(tissues == "Tissues")
tissues_quant_by_point_D14_blood$mouse <- as_factor(tissues_quant_by_point_D14_blood$mouse)
tissues_quant_by_point_D14_tissues$mouse <- as_factor(tissues_quant_by_point_D14_tissues$mouse)
wilcox.test(tissues_quant_by_point_D14_blood$total, tissues_quant_by_point_D14_tissues$total, p.adjust.methods = "BH")
