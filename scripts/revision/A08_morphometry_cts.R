# A8 — CTS occupancy / turnover vs morphometry and trophic state
source("scripts/revision/_bootstrap.R")
shared <- revision_bootstrap()

message("=== A8 morphometric contingency ===")

meta <- shared$lake_meta
trophic <- read.csv("data/revision/lake_trophic_state.csv", stringsAsFactors = FALSE) %>%
  dplyr::mutate(lake = revision_normalize_lake_names(lake))

# Prefer k=6 CTS1 occupancy if available; else compute quickly
cts_path <- "outputs/revision/A7_sample_cts6_assignments.csv"
if (!file.exists(cts_path)) {
  message("A7 assignments missing; running A7 first...")
  source("scripts/revision/A07_split_cts5_cts6.R")
}
cts <- readr::read_csv(cts_path, show_col_types = FALSE)

cts1_occ <- cts %>%
  dplyr::filter(age_ce > 0) %>%
  dplyr::group_by(lake) %>%
  dplyr::summarise(
    n_samples = dplyr::n(),
    n_CTS1 = sum(amd_clusts == "CTS1"),
    prop_CTS1 = mean(amd_clusts == "CTS1"),
    .groups = "drop"
  )

# Turnover magnitude from existing timing comparison if present
timing_path <- "outputs/lake_gam_change_timing_comparison.csv"
if (file.exists(timing_path)) {
  timing <- readr::read_csv(timing_path, show_col_types = FALSE)
} else if (file.exists("outputs/revision/A2_onset_comparison_by_resolution.csv")) {
  timing <- readr::read_csv("outputs/revision/A2_onset_comparison_by_resolution.csv", show_col_types = FALSE) %>%
    dplyr::filter(resolution == "species")
} else {
  stop("Need lake timing comparison CSV (run main pipeline or A2 first)")
}

lake_df <- timing %>%
  dplyr::transmute(
    lake,
    onset_leader,
    total_dca_change_Producers,
    total_dca_change_Consumers,
    magnitude_diff,
    mean_turnover = (total_dca_change_Producers + total_dca_change_Consumers) / 2
  ) %>%
  dplyr::left_join(cts1_occ, by = "lake") %>%
  dplyr::left_join(meta %>% dplyr::select(lake, island, alt, area, zmax), by = "lake") %>%
  dplyr::left_join(trophic, by = "lake")

revision_write_csv(lake_df, "outputs/revision/A8_lake_morpho_cts_turnover.csv")

cors <- purrr::map_dfr(
  c("prop_CTS1", "mean_turnover", "total_dca_change_Producers", "total_dca_change_Consumers"),
  function(resp) {
    purrr::map_dfr(c("zmax", "area", "alt", "tp_ug_l"), function(pred) {
      if (!pred %in% names(lake_df)) return(NULL)
      dplyr::bind_rows(
        revision_cor_pair(lake_df[[pred]], lake_df[[resp]], "pearson"),
        revision_cor_pair(lake_df[[pred]], lake_df[[resp]], "spearman")
      ) %>%
        dplyr::mutate(response = resp, predictor = pred)
    })
  }
)
revision_write_csv(cors, "outputs/revision/A8_morpho_correlations.csv")
print(cors)

# Deep-lake CTS1 note
deep_note <- lake_df %>%
  dplyr::arrange(dplyr::desc(zmax)) %>%
  dplyr::select(lake, zmax, prop_CTS1, mean_turnover, trophic_state)
revision_write_csv(deep_note, "outputs/revision/A8_depth_ranked_cts1.csv")

p <- ggplot(lake_df %>% dplyr::filter(is.finite(zmax)), aes(x = zmax, y = prop_CTS1, label = lake)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 20) +
  theme_minimal(base_size = 11) +
  labs(title = "A8: CTS1 occupancy vs Zmax", x = "Zmax (m)", y = "Proportion CTS1 samples")
# ggrepel optional
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  p <- ggplot(lake_df %>% dplyr::filter(is.finite(zmax)), aes(x = zmax, y = prop_CTS1)) +
    geom_point(size = 3) +
    geom_text(aes(label = lake), vjust = -0.7, size = 3) +
    theme_minimal(base_size = 11) +
    labs(title = "A8: CTS1 occupancy vs Zmax", x = "Zmax (m)", y = "Proportion CTS1 samples")
}
ggsave("outputs/revision/figures/A8_cts1_vs_zmax.png", p, width = 7, height = 5, dpi = 150)

p2 <- ggplot(lake_df %>% dplyr::filter(is.finite(zmax)), aes(x = zmax, y = mean_turnover, label = lake)) +
  geom_point(size = 3, aes(colour = onset_leader)) +
  geom_text(vjust = -0.7, size = 3) +
  theme_minimal(base_size = 11) +
  labs(title = "A8: Mean turnover vs Zmax", x = "Zmax (m)", y = "Mean |DCA change|")
ggsave("outputs/revision/figures/A8_turnover_vs_zmax.png", p2, width = 7.5, height = 5, dpi = 150)

message("A8 complete (fill Table S3 values in data/revision/lake_trophic_state.csv for TP tests)")
