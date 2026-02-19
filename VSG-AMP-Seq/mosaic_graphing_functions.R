# VSG-AMP-Seq functions
# Smith 2024
# Author: Jaclyn Smith

library(tidyverse)
library(qPCRr)

#negate the %in% function
"%notin%" <- Negate("%in%")

get_all_CO_events <- function(file_name, num_errors_allowed = 5, day_order, genotype_order, mouse_order,
                              primer_order, family = c("VSG-4156", "VSG-3110", "VSG-228", "VSG-2986", "VSG-7358")){
  indv_cross <-readr::read_csv(file_name,
                               col_names = c("mouse", "day", "genotype", "primer", "direction", "read_count",
                                             "donor_VSG", "read_start", "event_start", "event_end", "type_of_event", "UMI"),
                               na = c("", "NA"))
  indv_cross$genotype[is.na(indv_cross$genotype)] <- "WT"
  indv_cross$genotype <- stringr::str_replace(indv_cross$genotype, "muMt", "muMT")
  indv_cross$donor_VSG <- indv_cross$donor_VSG %>%
    stringr::str_replace("Tbb1125", "")
  
  weed_out_errors <- indv_cross %>%
    dplyr::group_by(read_count, donor_VSG) %>%
    dplyr::filter(type_of_event != "donor_to_target") %>%
    dplyr::filter(type_of_event != "target_to_donor") %>%
    dplyr::filter(type_of_event != "target_to_donor_edge") %>%
    dplyr::filter(type_of_event != "donor_to_target_edge") %>%
    dplyr::summarise(count = n()) %>%
    dplyr::filter(count > num_errors_allowed)
  
  trimmed_COs <- dplyr::anti_join(indv_cross, weed_out_errors, by = "read_count")
  
  # # remove ambiguous donors 
  trimmed_stat <- trimmed_COs %>%
    dplyr::group_by(read_count, donor_VSG) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(read_count) %>%
    dplyr::summarise(count_donors = n()) %>%
    dplyr::filter(count_donors > 1)
  
  # replace ambiguous recombination events with average donor position and remove errors
  # since the known donor is not clear errors associated with these events have been consolidated to ambiguous
  multiple_donors_removed <- anti_join(trimmed_COs, trimmed_stat, by = "read_count")
  ambiguous <- setdiff(trimmed_COs, multiple_donors_removed)%>%
    mutate(family_status = case_when(
      donor_VSG %in% family ~ "family",
      TRUE ~ "Ambiguous"
    ))
  
  # deal with ambiguous cross overs
  ambig_COs <- ambiguous %>%
    filter(type_of_event %in% c("target_to_donor", "donor_to_target", "target_to_donor_edge", "donor_to_target_edge"))
  
  actually_unambig_COs <- ambig_COs %>%
    group_by(read_count) %>%
    summarise(count = n()) %>%
    filter(count == 1)
  View(actually_unambig_COs)
  
  # when there are family and non-family identified donors, all are grouped under fully ambiguous
  ambig_COs_summary <- ambig_COs %>%
    group_by(read_count, family_status) %>%
    summarise(average_needed = n()) %>%
    ungroup() %>%
    group_by(read_count) %>%
    summarise(both_family = n()) %>%
    filter(both_family > 1)
  
  if (length(ambig_COs_summary$read_count) > 0){
    ambig_COs <- ambig_COs %>%
      mutate(family_status = case_when(
        read_count %in% ambig_COs_summary$read_count ~ "Ambiguous",
        read_count %in% actually_unambig_COs$read_count ~ donor_VSG,
        TRUE ~ family_status))
  }
  
  # averge_all_recomb_identified
  ambig_COs <- ambig_COs %>%
    group_by(mouse, day, genotype, primer, direction, read_start, type_of_event, read_count, UMI) %>%
    summarise(event_start = mean(event_start), event_end = mean(event_end), donor_VSG = family_status) %>%
    ungroup() %>%
    group_by(mouse, day, genotype, primer, direction, read_start, event_start, event_end, donor_VSG, type_of_event,
             read_count, UMI) %>%
    summarise(ambigcount = n()) %>%
    select(-ambigcount)
  
  # View(ambig_COs)
  # deal with ambiguous errors
  ambig_errors <- ambiguous %>%
    filter(type_of_event %in% c("random_error", "random_target", "random_donor"))
  
  unambig_donor <- ambiguous %>%
    filter(type_of_event %in% c("target_to_donor", "donor_to_target", "target_to_donor_edge", "donor_to_target_edge")) %>%
    filter(read_count %in% actually_unambig_COs$read_count)
  
  change_status_to_ambiguous <- ambig_errors %>%
    group_by(read_count, family_status) %>%
    summarise(count = n()) %>%
    ungroup() %>%
    group_by(read_count) %>%
    summarise(count = n()) %>%
    filter(count > 1)
  View(change_status_to_ambiguous)

  # remove errors when they do not match the unambiguous donor
  # update family to ambiguous if multiple donors identified
  if (length(actually_unambig_COs$read_count) > 0){
    ambig_errors_unambig <- ambig_errors %>%
      mutate(family_status = case_when(
        read_count %in% actually_unambig_COs$read_count & donor_VSG == unambig_donor$donor_VSG ~ donor_VSG,
        read_count %in% actually_unambig_COs$read_count ~ "remove",
        read_count %in% change_status_to_ambiguous$read_count ~"Ambiguous",
        TRUE ~ family_status
      )) %>%
      filter(family_status != "remove")
    
  }
  else{
    ambig_errors_unambig <- ambig_errors %>%
      mutate(family_status = case_when(
        read_count %in% change_status_to_ambiguous$read_count ~"Ambiguous",
        TRUE ~ family_status
      )) %>%
      filter(family_status != "remove")
  }
  ambig_associated_donors <- ambig_errors_unambig %>%
    group_by(mouse, day, genotype, primer, direction, read_count, donor_VSG) %>%
    summarise(count = n()) %>%
    ungroup() %>%
    group_by(read_count) %>%
    summarise(VSGcount = n())
  
  ambig_associated_errors <- ambig_errors_unambig %>%
    group_by(mouse, day, genotype, primer, direction, read_count,
             type_of_event, read_start, event_start, event_end, UMI)  %>%
    summarise(count = n(), donor_VSG = family_status)
  
  ambig_errors_sum <- left_join(ambig_associated_errors, ambig_associated_donors) %>%
    filter(count == VSGcount) %>%
    select(-count, -VSGcount) %>%
    group_by(mouse, day, genotype, primer, direction, read_count,
             type_of_event, read_start, event_start, event_end, UMI, donor_VSG) %>%
    summarise(count = n()) %>%
    select(-count)
  
  final_trim <- bind_rows(multiple_donors_removed, ambig_COs, ambig_errors_sum)
  
  
  updated_singleCOs <- final_trim %>%
    dplyr::mutate(target_start = read_start + event_start) %>%
    dplyr::mutate(target_end = read_start + event_end) %>%
    dplyr::mutate(avg_pos = (target_start + target_end)/2)
  
  updated_singleCOs <- qPCRr::reorder_samples(updated_singleCOs, "day", day_order)
  updated_singleCOs <- qPCRr::reorder_samples(updated_singleCOs, "mouse", mouse_order)
  select_genotypes = intersect(genotype_order, updated_singleCOs$genotype)
  updated_singleCOs <- qPCRr::reorder_samples(updated_singleCOs, "genotype", select_genotypes)
  select_primers = intersect(primer_order, updated_singleCOs$primer)
  updated_singleCOs <- qPCRr::reorder_samples(updated_singleCOs, "primer", select_primers)

  
  return(updated_singleCOs)
}

color_VSG_vector = c("VSG-228" = "turquoise3", "VSG-2986" = "slateblue4", "VSG-3110" = "steelblue3", "VSG-7358" = "darkorchid1",
                     "VSG-4156" = "darkslategray4", "VSG-1128" = "lightpink2", "VSG-479" = "indianred4", "other VSGs" = "grey80", "family" = "dodgerblue3")

graph_all_events <- function(single_events, color_var = "day", wrap_by = "mouse", xlab = "Average Position",
                             ylab = "Count", leg = "Day", D14 = FALSE, select_VSGs = FALSE, remove_VSGs = FALSE,
                             vsg_interest_list = c(), vsg_remove_list = c(), leg_pos = "right", multiple = FALSE, merge = FALSE, merge_columns = c(),
                             merge_order = c(), add_cutline = FALSE, bin_width = 25, cut_line = FALSE, cut_line_names = c(),
                             cut_line_locations = c(), cut_line_colors = c(), adjust_xaxis = FALSE, color_scheme = TRUE, vert_num = 6,
                             log_scale = FALSE, summary_graph = FALSE, tissues = FALSE, by_mouse_only = FALSE, mouse_plot = FALSE, rep_pic = FALSE,
                             seq_type = FALSE){
  if (summary_graph) {
    all_single_events <- single_events %>%
      dplyr::group_by(donor_VSG, avg_pos, type_of_event) %>%
      dplyr::summarise(count = n()) 
  } else if (by_mouse_only){
    all_single_events <- single_events %>%
      mutate(tissues = case_when(
        genotype == "Blood" ~ "Blood",
        TRUE ~ "Tissues")) %>%
      dplyr::group_by(mouse, day, donor_VSG, avg_pos, type_of_event, tissues) %>%
      dplyr::summarise(count = n()) 
  } else if (rep_pic){
    all_single_events <- single_events %>%
      dplyr::group_by(genotype, donor_VSG, avg_pos, type_of_event) %>%
      dplyr::summarise(count = n()) 
  } else if (seq_type) {
    all_single_events <- single_events %>%
      dplyr::group_by(mouse, day, genotype, donor_VSG, avg_pos, type_of_event, seq_type) %>%
      dplyr::summarise(count = n())
  } else {
    all_single_events <- single_events %>%
      dplyr::group_by(mouse, day, genotype, donor_VSG, avg_pos, type_of_event) %>%
      dplyr::summarise(count = n()) 
  }
  
  
  if (multiple) {
    all_single_events <- single_events %>%
      dplyr::group_by(mouse, day, genotype, read_count, type_of_event, donor_VSG, avg_pos) %>%
      dplyr::summarise(count = n())
  }
  
  if (D14) {
    all_single_events <- all_single_events %>%
      dplyr::filter(day == "D14")
  }
  
  if (select_VSGs) {
    all_single_events <- all_single_events %>%
      dplyr::filter(donor_VSG %in% vsg_interest_list)
    all_single_events <- qPCRr::reorder_samples(all_single_events, "donor_VSG", vsg_interest_list)
  }
  
  if (remove_VSGs) {
    all_single_events <- all_single_events %>%
      dplyr::filter(donor_VSG %notin% vsg_remove_list)
  }
  
  if (merge) {
    new_column = paste(merge_columns, collapse = "_")
    all_single_events <- all_single_events %>%
      tidyr::unite(!!new_column, merge_columns)
    all_single_events <- qPCRr::reorder_samples(all_single_events, new_column, merge_order)
  }
  
  if (tissues) {
    all_single_events <- all_single_events %>%
      mutate(tissues = case_when(
        genotype == "Blood" ~ "Blood",
        TRUE ~ "Tissues"
      ))
  }
  
  View(all_single_events)
  
  graph_plot <- ggplot2::ggplot(all_single_events, aes(avg_pos)) +
    ggplot2::geom_histogram(ggplot2::aes_string(fill = color_var), binwidth = bin_width)
  
  if (mouse_plot) {
    graph_plot <- graph_plot + ggplot2::facet_wrap(wrap_by, nrow = vert_num, scales = "free_y")
  }
  else if (!summary_graph) {
    graph_plot <- graph_plot + ggplot2::facet_wrap(wrap_by, nrow = vert_num)
  }
  graph_plot <- graph_plot + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 20,
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
                                            legend.position = leg_pos,
                                            strip.text.x = ggplot2::element_text(size = 20),
                                            legend.title = ggplot2::element_text(size = 15,
                                                                                 face = "bold",
                                                                                 hjust = 0.5),
                                            legend.key = ggplot2::element_rect(fill = "white"),
                                            legend.text = ggplot2::element_text(size = 15)) +
    ggplot2::labs(x = xlab, y = ylab, fill = leg)
  
  
  if (cut_line) {
    added_lines <- data.frame(genotype = cut_line_names,
                              cut_site = cut_line_locations)
    added_lines <- qPCRr::reorder_samples(as_tibble(added_lines), "genotype", cut_line_names)
    graph_plot <- graph_plot +
      ggplot2::geom_vline(aes(xintercept = cut_site), data = added_lines, color = cut_line_colors)
  }
  
  if (adjust_xaxis) {
    graph_plot <- graph_plot +
      scale_x_continuous(breaks = seq(0,1600, by=400), limits=c(0, 1687))
  }
  
  if (log_scale) {
    graph_plot <- graph_plot +
      scale_y_continuous(trans = scales::log1p_trans())
  }
  
  if (color_scheme) {
    graph_plot <- graph_plot +
      scale_fill_manual(values = color_VSG_vector,
                        limits = c("VSG-228", "VSG-2986", "VSG-3110", "VSG-7358",
                                   "VSG-4156", "VSG-1128", "VSG-479", "other VSGs", "family"),
                        na.value = "grey80")
    
  }
  
  return(graph_plot)
}


