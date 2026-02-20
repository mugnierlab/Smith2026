# Supplemental Figure 9B
# Mouse Infection Extended Data
# Smith 2025

# experiment 1 cleaning data
raw_expt1 <- readxl::read_xlsx("/Users/jaclynsmith/labNotebook/monicaMugnierLabNotebook/JAC21-4/JAC21-4mouse_collectionExpt1.xlsx",
                               sheet = "Parasitemia")

expt1 <- raw_expt1 %>%
  dplyr::distinct() %>%
  purrr::set_names(raw_expt1[1,]) %>%
  dplyr::select(contains("D")) %>%
  dplyr:: slice(3:14)

multiplier <- expt1 %>%
  dplyr::slice(1, 7) %>%
  dplyr::mutate(genotype = c("WT", "muMT")) %>%
  tidyr::pivot_longer(!genotype, names_to = "day", values_to = "multiplier")
multiplier$multiplier <- as.numeric(multiplier$multiplier)


dates_expt1 <- expt1 %>%
  dplyr::slice(-1, -7) %>%
  dplyr::mutate(mouse = c("E1M6", "E1M7", "E1M8", "E1M9", "E1M10", "E1M1", "E1M2", "E1M3", "E1M4", "E1M5"),
                genotype = c("WT", "WT", "WT", "WT", "WT", "muMT", "muMT", "muMT", "muMT", "muMT"),
                experiment = "E1") %>%
  dplyr::na_if("NA") %>%
  tidyr::pivot_longer(!c(mouse, genotype, experiment), names_to = "day", values_to = "count") %>%
  tidyr::drop_na() %>%
  tidyr::separate(count, ",", into = c("one", "two", "three", "four"), fill = "right")

dates_expt1$one <- as.numeric(dates_expt1$one)
dates_expt1$two <- as.numeric(dates_expt1$two)
dates_expt1$three <- as.numeric(dates_expt1$three)
dates_expt1$four <- as.numeric(dates_expt1$four)

dates_expt1$average <- rowMeans(dates_expt1[,5:8], na.rm = TRUE)

expt1_final <- dates_expt1 %>%
  dplyr::select(mouse, genotype, day, experiment, average) %>%
  dplyr::full_join(multiplier, by = c("day", "genotype")) %>%
  dplyr::mutate(parasitemia = average*multiplier*(10^4)/(5/1000)) %>%
  dplyr::mutate(parasitemia = dplyr::case_when(
    parasitemia == 0 ~ 10000,
    TRUE ~ parasitemia
  )) %>%
  dplyr::mutate(parasitemia_log = log(parasitemia,10))

expt1_final$day <- expt1_final$day %>%
  stringr::str_replace("D", "")
expt1_final$day <- as.numeric(expt1_final$day)


###### clean expt 2

raw_expt2 <- readxl::read_xlsx("/Users/jaclynsmith/labNotebook/monicaMugnierLabNotebook/JAC21-4/JAC21-4mouse_collectionExpt2.xlsx",
                               sheet = "Parasitemia")

expt2 <- raw_expt2 %>%
  dplyr::distinct() %>%
  purrr::set_names(raw_expt1[1,]) %>%
  dplyr::select(contains("D")) %>%
  dplyr:: slice(3:14)

multiplier <- expt2 %>%
  dplyr::slice(1, 7) %>%
  dplyr::mutate(genotype = c("WT", "muMT")) %>%
  tidyr::pivot_longer(!genotype, names_to = "day", values_to = "multiplier")
multiplier$multiplier <- as.numeric(multiplier$multiplier)


dates_expt2 <- expt2 %>%
  dplyr::slice(-1, -7) %>%
  dplyr::mutate(mouse = c("E2M1", "E2M2", "E2M3", "E2M4", "E2M5", "E2M6", "E2M7", "E2M8", "E2M9", "E2M10"),
                genotype = c("WT", "WT", "WT", "WT", "WT", "muMT", "muMT", "muMT", "muMT", "muMT"),
                experiment = "E2") %>%
  dplyr::na_if("NA") %>%
  tidyr::pivot_longer(!c(mouse, genotype, experiment), names_to = "day", values_to = "count") %>%
  tidyr::drop_na() %>%
  tidyr::separate(count, ",", into = c("one", "two"), fill = "right")

dates_expt2$one <- as.numeric(dates_expt2$one)
dates_expt2$two <- as.numeric(dates_expt2$two)

dates_expt2$average <- rowMeans(dates_expt2[,5:6], na.rm = TRUE)

expt2_final <- dates_expt2 %>%
  dplyr::select(mouse, genotype, day, experiment, average) %>%
  dplyr::full_join(multiplier, by = c("day", "genotype")) %>%
  dplyr::mutate(parasitemia = average*multiplier*(10^4)/(5/1000)) %>%
  dplyr::mutate(parasitemia = dplyr::case_when(
    parasitemia == 0 ~ 10000,
    TRUE ~ parasitemia
  )) %>%
  dplyr::mutate(parasitemia_log = log(parasitemia,10))

expt2_final$day <- expt2_final$day %>%
  stringr::str_replace("D", "")
expt2_final$day <- as.numeric(expt2_final$day)

#####EXPT3 CLEANING

raw_expt3 <- readxl::read_xlsx("/Users/jaclynsmith/labNotebook/monicaMugnierLabNotebook/JAC21-4/JAC21-4mouse_collectionExpt3.xlsx",
                               sheet = "Parasitemia")

expt3 <- raw_expt3 %>%
  dplyr::distinct() %>%
  purrr::set_names(raw_expt3[1,]) %>%
  dplyr::select(contains("D")) %>%
  dplyr:: slice(3:14)

multiplier <- expt3 %>%
  dplyr::slice(1, 7) %>%
  dplyr::mutate(genotype = c("WT", "muMT")) %>%
  tidyr::pivot_longer(!genotype, names_to = "day", values_to = "multiplier")
multiplier$multiplier <- as.numeric(multiplier$multiplier)


dates_expt3 <- expt3 %>%
  dplyr::slice(-1, -7) %>%
  dplyr::mutate(mouse = c("E3M6", "E3M7", "E3M8", "E3M9", "E3M10", "E3M1", "E3M2", "E3M3", "E3M4", "E3M5"),
                genotype = c("WT", "WT", "WT", "WT", "WT", "muMT", "muMT", "muMT", "muMT", "muMT"),
                experiment = "E3") %>%
  dplyr::na_if("NA") %>%
  tidyr::pivot_longer(!c(mouse, genotype, experiment), names_to = "day", values_to = "count") %>%
  tidyr::drop_na() %>%
  tidyr::separate(count, ",", into = c("one", "two", "three"), fill = "right")

dates_expt3$one <- as.numeric(dates_expt3$one)
dates_expt3$two <- as.numeric(dates_expt3$two)
dates_expt3$three <- as.numeric(dates_expt3$three)


dates_expt3$average <- rowMeans(dates_expt3[,5:7], na.rm = TRUE)

expt3_final <- dates_expt3 %>%
  dplyr::select(mouse, genotype, day, experiment, average) %>%
  dplyr::full_join(multiplier, by = c("day", "genotype")) %>%
  dplyr::mutate(parasitemia = average*multiplier*(10^4)/(5/1000)) %>%
  dplyr::mutate(parasitemia = dplyr::case_when(
    parasitemia == 0 ~ 10000,
    TRUE ~ parasitemia
  )) %>%
  dplyr::mutate(parasitemia_log = log(parasitemia,10))

expt3_final$day <- expt3_final$day %>%
  stringr::str_replace("D", "")
expt3_final$day <- as.numeric(expt3_final$day)



############################

all_parasitemia <- bind_rows(expt1_final, expt2_final, expt3_final)




total_parasitemia_graph <- graph_parasitemia(all_parasitemia, genotype_order = c("WT", "muMT"), remove_mice = TRUE,
                                             mouse_interest_list = c("E1M9", "E1M6",  "E1M1", "E1M2", "E1M3", "E2M10", "E2M4", "E2M5",
                                                                     "E3M1", "E2M3",
                                                                     "E3M2", "E3M5", "E3M3", "EMM1", "EMM2", "EMM3", "EMM4", "EMM5", "EMM6",
                                                                     "EMM7", "EMM8", "EMM9", "EMM10", "E2M4",
                                                                     "E2M5", "E3M8", "E3M9", "E3M10", "E1M6"), wrap_by = "experiment", color_var = "genotype") +
  scale_color_manual(values = c("darkslategray", "darksalmon"), limits = c("WT", "muMT"), na.value = "grey80")
