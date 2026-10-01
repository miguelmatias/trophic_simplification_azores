# A9 — Lake onset order vs fish-stocking dates
source("scripts/revision/_bootstrap.R")
shared <- revision_bootstrap()

message("=== A9 onset vs stocking ===")

stock <- read.csv("data/revision/fish_stocking.csv", stringsAsFactors = FALSE, na.strings = c("NA", "")) %>%
  dplyr::mutate(
    lake = revision_normalize_lake_names(lake),
    first_stocking_ce = as.numeric(first_stocking_ce)
  )

timing_path <- "outputs/lake_gam_change_timing_comparison.csv"
if (!file.exists(timing_path)) {
  a2_path <- "outputs/revision/A2_onset_comparison_by_resolution.csv"
  if (file.exists(a2_path)) {
    timing <- readr::read_csv(a2_path, show_col_types = FALSE) %>%
      dplyr::filter(resolution == "species")
  } else {
    stop("Missing timing comparison; run main script S7 chunk or A2 first")
  }
} else {
  timing <- readr::read_csv(timing_path, show_col_types = FALSE)
}

out <- timing %>%
  dplyr::left_join(stock, by = "lake") %>%
  dplyr::mutate(
    consumer_onset = first_onset_after_min_ce_Consumers,
    producer_onset = first_onset_after_min_ce_Producers,
    consumer_after_stocking = dplyr::case_when(
      is.na(first_stocking_ce) | is.na(consumer_onset) ~ NA,
      consumer_onset >= first_stocking_ce ~ TRUE,
      TRUE ~ FALSE
    ),
    producer_after_stocking = dplyr::case_when(
      is.na(first_stocking_ce) | is.na(producer_onset) ~ NA,
      producer_onset >= first_stocking_ce ~ TRUE,
      TRUE ~ FALSE
    ),
    lag_consumer_minus_stocking = consumer_onset - first_stocking_ce,
    lag_producer_minus_stocking = producer_onset - first_stocking_ce
  )

revision_write_csv(out, "outputs/revision/A9_onset_vs_stocking.csv")

summary_tbl <- out %>%
  dplyr::summarise(
    n_lakes = dplyr::n(),
    n_with_stocking_date = sum(is.finite(first_stocking_ce)),
    n_producers_earlier = sum(onset_leader == "Producers earlier", na.rm = TRUE),
    n_consumers_earlier = sum(onset_leader == "Consumers earlier", na.rm = TRUE),
    n_similar = sum(onset_leader == "Similar timing", na.rm = TRUE),
    consumer_first_lakes = paste(lake[onset_leader == "Consumers earlier"], collapse = "; ")
  )
revision_write_csv(summary_tbl, "outputs/revision/A9_onset_stocking_summary.csv")
print(out %>% dplyr::select(lake, onset_leader, producer_onset, consumer_onset, first_stocking_ce, lag_consumer_minus_stocking))

p <- ggplot(out, aes(x = reorder(lake, producer_onset))) +
  geom_point(aes(y = producer_onset, colour = "Producers"), size = 3) +
  geom_point(aes(y = consumer_onset, colour = "Consumers"), size = 3) +
  geom_point(
    data = out %>% dplyr::filter(is.finite(first_stocking_ce)),
    aes(y = first_stocking_ce, shape = "Stocking"),
    size = 3, colour = "black"
  ) +
  coord_flip() +
  theme_minimal(base_size = 11) +
  labs(
    title = "A9: Onset order vs fish stocking",
    y = "Age (CE)", x = NULL, colour = "Onset", shape = NULL
  )
ggsave("outputs/revision/figures/A9_onset_vs_stocking.png", p, width = 8, height = 5, dpi = 150)

message("A9 complete (add stocking dates for non-Azul lakes in data/revision/fish_stocking.csv)")
