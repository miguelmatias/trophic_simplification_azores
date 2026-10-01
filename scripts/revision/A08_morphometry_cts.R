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

# Ordinal trophic state for Spearman (oligo < meso < eu); NA if unknown
trophic_ord_levels <- c(
  "oligotrophic", "oligo-mesotrophic", "mesotrophic",
  "meso-eutrophic", "eutrophic", "hypereutrophic"
)

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
  dplyr::left_join(trophic, by = "lake") %>%
  dplyr::mutate(
    # Prata is now a peatland (Zmax = 0 / unavailable) — exclude from depth tests
    zmax_for_depth = dplyr::if_else(
      lake == "Prata" | !is.finite(zmax) | zmax <= 0,
      NA_real_,
      zmax
    ),
    trophic_ord = as.integer(factor(
      tolower(trimws(as.character(trophic_state))),
      levels = trophic_ord_levels,
      ordered = TRUE
    ))
  )

revision_write_csv(lake_df, "outputs/revision/A8_lake_morpho_cts_turnover.csv")

# Depth / morphometry / trophic correlations (n = 8 lakes with finite Zmax > 0)
cors <- purrr::map_dfr(
  c("prop_CTS1", "mean_turnover", "total_dca_change_Producers", "total_dca_change_Consumers"),
  function(resp) {
    purrr::map_dfr(
      c("zmax_for_depth", "area", "alt", "tp_ug_l", "trophic_ord"),
      function(pred) {
        if (!pred %in% names(lake_df)) return(NULL)
        dplyr::bind_rows(
          revision_cor_pair(lake_df[[pred]], lake_df[[resp]], "pearson"),
          revision_cor_pair(lake_df[[pred]], lake_df[[resp]], "spearman")
        ) %>%
          dplyr::mutate(
            response = resp,
            predictor = dplyr::recode(pred, zmax_for_depth = "zmax", trophic_ord = "trophic_state")
          )
      }
    )
  }
)
revision_write_csv(cors, "outputs/revision/A8_morpho_correlations.csv")
print(cors)

# Focused Spearman table for reviewer response (CTS1 + turnover × Zmax/area/trophic)
focus <- cors %>%
  dplyr::filter(
    method == "spearman",
    response %in% c("prop_CTS1", "mean_turnover"),
    predictor %in% c("zmax", "area", "trophic_state", "tp_ug_l")
  ) %>%
  dplyr::arrange(response, predictor)
revision_write_csv(focus, "outputs/revision/A8_spearman_focus.csv")
print(focus)

# Deep-lake CTS1 note (Prata excluded from depth ranking)
deep_note <- lake_df %>%
  dplyr::filter(is.finite(zmax_for_depth)) %>%
  dplyr::arrange(dplyr::desc(zmax_for_depth)) %>%
  dplyr::select(lake, zmax, prop_CTS1, mean_turnover, trophic_state)
revision_write_csv(deep_note, "outputs/revision/A8_depth_ranked_cts1.csv")

prata_note <- tibble::tibble(
  lake = "Prata",
  reason = "Excluded from Zmax analyses: now a peatland (Zmax = 0 / unavailable in metadata).",
  n_morpho_lakes = sum(is.finite(lake_df$zmax_for_depth))
)
revision_write_csv(prata_note, "outputs/revision/A8_prata_exclusion_note.csv")

plot_df <- lake_df %>% dplyr::filter(is.finite(zmax_for_depth))

rho_cts1_z <- focus %>%
  dplyr::filter(response == "prop_CTS1", predictor == "zmax") %>%
  dplyr::slice(1)
rho_turn_z <- focus %>%
  dplyr::filter(response == "mean_turnover", predictor == "zmax") %>%
  dplyr::slice(1)

fmt_rho <- function(r) {
  if (nrow(r) == 0 || is.na(r$estimate[1])) return("Spearman ρ = NA")
  sprintf("Spearman ρ = %.2f, p = %.3f (n = %d)", r$estimate[1], r$p.value[1], r$n[1])
}

p <- ggplot(plot_df, aes(x = zmax_for_depth, y = prop_CTS1, label = lake)) +
  geom_point(size = 3) +
  theme_minimal(base_size = 11) +
  labs(
    title = "A8: CTS1 occupancy vs Zmax",
    subtitle = paste0(
      "CTS1 = euplanctonic-dominated (k=6). ", fmt_rho(rho_cts1_z),
      ". Prata excluded (peatland, Zmax = 0)."
    ),
    x = "Zmax (m)",
    y = "Proportion CTS1 samples"
  )
if (requireNamespace("ggrepel", quietly = TRUE)) {
  p <- p + ggrepel::geom_text_repel(size = 3, max.overlaps = 20)
} else {
  p <- p + geom_text(vjust = -0.7, size = 3)
}
ggsave("outputs/revision/figures/A8_cts1_vs_zmax.png", p, width = 7.5, height = 5.2, dpi = 150)

p2 <- ggplot(plot_df, aes(x = zmax_for_depth, y = mean_turnover, label = lake)) +
  geom_point(size = 3, aes(colour = onset_leader)) +
  geom_text(vjust = -0.7, size = 3) +
  theme_minimal(base_size = 11) +
  labs(
    title = "A8: Mean turnover vs Zmax",
    subtitle = paste0(
      fmt_rho(rho_turn_z),
      ". Funda is a high-turnover outlier; Prata excluded (peatland)."
    ),
    x = "Zmax (m)",
    y = "Mean |DCA change|"
  )
ggsave("outputs/revision/figures/A8_turnover_vs_zmax.png", p2, width = 7.5, height = 5.2, dpi = 150)

message("A8 complete")
message("Prata excluded from Zmax tests (peatland / Zmax unavailable).")
message("Fill Table S3 values in data/revision/lake_trophic_state.csv for trophic/TP tests with n≥4.")
