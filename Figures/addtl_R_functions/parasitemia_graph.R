graph_parasitemia <- function(parasitemia, color_var = "genotype", wrap_by = "experiment", xlab = "Day",
                              ylab = "Parasitemia (parasites/mL)", leg = "Genotype", wrap = TRUE, mouse_order,
                              genotype_order = c("WT", "muMT"), remove_mice = FALSE, select_mice = FALSE, mouse_interest_list = c()) {
  
  if (remove_mice) {
    parasitemia <- parasitemia %>%
      dplyr::filter(mouse %notin% mouse_interest_list)
  }

  
  
  parasitemia <- qPCRr::reorder_samples(parasitemia, "genotype", genotype_order)

  parasite_graph <- ggplot(parasitemia, aes_string(x = "day", y = "parasitemia_log", group = "genotyper", color = color_var)) + 
    geom_jitter(width = 0.5) +
    geom_line() +
    facet_wrap(wrap_by) +
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
    labs(x = xlab, y = ylab, legend = leg) +
    scale_x_continuous(breaks = seq(4, 22, by = 2)) 
  return(parasite_graph)
  
}
