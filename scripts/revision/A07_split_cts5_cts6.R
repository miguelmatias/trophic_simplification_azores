# A7 — Compare k=5 (independent AMD), k=6, and manuscript k=6→k=5 merge
source("scripts/revision/_bootstrap.R")
shared <- revision_bootstrap()

message("=== A7 CTS schemes: independent k=5, k=6, manuscript merge ===")

# Rebuild normalized abundances from df_fgroups (same as main_script)
fg <- shared$df_fgroups %>%
  dplyr::filter(is.finite(age_ce)) %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(revision_guild_cols), ~ replace(.x, is.na(.x), 0)))

norm_abund <- dplyr::bind_cols(
  fg %>% dplyr::select(lake, core_depth_id, age_ce),
  {
    prod <- as.matrix(fg[c("high_profile", "low_profile", "motile", "euplanktonic")])
    cons <- as.matrix(fg[c("algivore", "detritivore", "plantivore", "predator")])
    prod <- prod / pmax(rowSums(prod), 1e-8)
    cons <- cons / pmax(rowSums(cons), 1e-8)
    as.data.frame(cbind(prod, cons))
  }
) %>%
  replace(is.na(.), 0) %>%
  dplyr::ungroup()

# Match main_script AMD iterations at fixed k (asymptote there typically yields k = 6)
amd_iterations <- 5000L
set.seed(42)
message("Clustering k=6 (single partition for k=6 + manuscript merge)...")
assign_k6 <- revision_run_amd_raw(norm_abund, k = 6L, iterations = amd_iterations)
cts6 <- revision_relabel_cts_euplanctonic(norm_abund, assign_k6, k = 6L)

merged <- revision_manuscript_k6_merged_to_k5(norm_abund, assign_k6)
cts5_merged <- merged$data

message("Clustering k=5 (independent AMD)...")
assign_k5 <- revision_run_amd_raw(norm_abund, k = 5L, iterations = amd_iterations)
cts5_indep <- revision_relabel_cts_euplanctonic(norm_abund, assign_k5, k = 5L)

scheme_labels <- c(
  k5_independent_amd = "k=5 (independent AMD)",
  k6 = "k=6",
  k5_merged_manuscript = "k=5 (merged from k=6, manuscript)"
)

revision_write_csv(merged$merge_map, "outputs/revision/A7_k6_manuscript_merge_map.csv")

# Ranking keys (euplanctonic-ordered) for audit
rank_key <- function(df, label) {
  df %>%
    dplyr::group_by(amd_clusts) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean_euplanktonic = mean(euplanktonic, na.rm = TRUE),
      mean_consumer = mean(algivore + detritivore + plantivore + predator, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(amd_clusts) %>%
    dplyr::mutate(scheme = label)
}

revision_write_csv(
  dplyr::bind_rows(
    rank_key(cts5_indep, scheme_labels[["k5_independent_amd"]]),
    rank_key(cts6, scheme_labels[["k6"]]),
    rank_key(cts5_merged, scheme_labels[["k5_merged_manuscript"]])
  ),
  "outputs/revision/A7_cts_euplanctonic_rank_key.csv"
)

profiles <- function(df, label) {
  df %>%
    dplyr::group_by(amd_clusts) %>%
    dplyr::summarise(
      n = dplyr::n(),
      dplyr::across(dplyr::all_of(revision_guild_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(scheme = label)
}

prof <- dplyr::bind_rows(
  profiles(cts5_indep, scheme_labels[["k5_independent_amd"]]),
  profiles(cts6, scheme_labels[["k6"]]),
  profiles(cts5_merged, scheme_labels[["k5_merged_manuscript"]])
)
revision_write_csv(prof, "outputs/revision/A7_cts_guild_profiles.csv")

div_sum <- dplyr::bind_rows(
  revision_cts_diversity_summary(cts5_indep, scheme_label = scheme_labels[["k5_independent_amd"]]),
  revision_cts_diversity_summary(cts6, scheme_label = scheme_labels[["k6"]]),
  revision_cts_diversity_summary(cts5_merged, scheme_label = scheme_labels[["k5_merged_manuscript"]])
)
revision_write_csv(div_sum, "outputs/revision/A7_cts_diversity_by_scheme.csv")

occ <- dplyr::bind_rows(
  revision_cts_occupancy(cts5_merged) %>% dplyr::mutate(scheme = scheme_labels[["k5_merged_manuscript"]]),
  revision_cts_occupancy(cts5_indep) %>% dplyr::mutate(scheme = scheme_labels[["k5_independent_amd"]]),
  revision_cts_occupancy(cts6) %>% dplyr::mutate(scheme = scheme_labels[["k6"]])
)
revision_write_csv(occ, "outputs/revision/A7_cts_occupancy_by_lake.csv")

temporal <- function(df, label, k) {
  levels_k <- paste0("CTS", seq_len(k))
  df %>%
    dplyr::filter(age_ce > 0) %>%
    dplyr::mutate(
      age_bin = floor(age_ce / 30) * 30,
      amd_clusts = factor(as.character(amd_clusts), levels = levels_k)
    ) %>%
    dplyr::count(age_bin, amd_clusts, name = "n", .drop = FALSE) %>%
    dplyr::group_by(age_bin) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(scheme = label) %>%
    dplyr::rename(age_ce = age_bin)
}

temp <- dplyr::bind_rows(
  temporal(cts5_indep, scheme_labels[["k5_independent_amd"]], 5L),
  temporal(cts6, scheme_labels[["k6"]], 6L),
  temporal(cts5_merged, scheme_labels[["k5_merged_manuscript"]], 5L)
)
bin_sums <- temp %>%
  dplyr::group_by(scheme, age_ce) %>%
  dplyr::summarise(prop_sum = sum(prop), .groups = "drop")
stopifnot(all(abs(bin_sums$prop_sum - 1) < 1e-9))
revision_write_csv(temp, "outputs/revision/A7_cts_temporal_proportions.csv")
revision_write_csv(bin_sums, "outputs/revision/A7_cts_temporal_bin_sums.csv")

pick_cts <- function(df, scheme) {
  df %>%
    dplyr::select(lake, core_depth_id, age_ce, amd_clusts) %>%
    dplyr::mutate(scheme = scheme)
}

early_stability <- dplyr::bind_rows(
  pick_cts(cts5_indep, scheme_labels[["k5_independent_amd"]]),
  pick_cts(cts6, scheme_labels[["k6"]]),
  pick_cts(cts5_merged, scheme_labels[["k5_merged_manuscript"]])
) %>%
  dplyr::filter(age_ce > 0, age_ce < 750) %>%
  dplyr::count(scheme, amd_clusts) %>%
  dplyr::group_by(scheme) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup()
revision_write_csv(early_stability, "outputs/revision/A7_early_phase_cts_composition.csv")
print(early_stability)

# Sample assignments (manuscript k=5 = merged; A1x / Fig. 3e)
revision_write_csv(
  cts6 %>% dplyr::select(lake, core_depth_id, age_ce, amd_clusts),
  "outputs/revision/A7_sample_cts6_assignments.csv"
)
revision_write_csv(
  cts5_merged %>% dplyr::select(lake, core_depth_id, age_ce, amd_clusts),
  "outputs/revision/A7_sample_cts5_assignments.csv"
)
revision_write_csv(
  cts5_merged %>% dplyr::select(lake, core_depth_id, age_ce, amd_clusts),
  "outputs/revision/A7_sample_cts5_merged_assignments.csv"
)
revision_write_csv(
  cts5_indep %>% dplyr::select(lake, core_depth_id, age_ce, amd_clusts),
  "outputs/revision/A7_sample_cts5_independent_amd_assignments.csv"
)

cts_cols_6 <- revision_cts_colours(6L)
cts_cols_5 <- revision_cts_colours(5L)

p <- ggplot(temp, aes(x = age_ce, y = prop, fill = amd_clusts)) +
  geom_area(position = "stack", alpha = 0.9) +
  scale_fill_manual(values = cts_cols_6, name = "amd_clusts", drop = FALSE) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  facet_wrap(~scheme, ncol = 1) +
  theme_minimal(base_size = 11) +
  labs(
    title = "A7: CTS temporal occupancy (three clustering schemes)",
    x = "Age (CE)",
    y = "Proportion"
  )
ggsave("outputs/revision/figures/A7_cts_temporal.png", p, width = 9, height = 10, dpi = 150)

prof_long <- prof %>%
  tidyr::pivot_longer(dplyr::all_of(revision_guild_cols), names_to = "guild", values_to = "mean_rel") %>%
  dplyr::mutate(
    scheme = factor(scheme, levels = unname(scheme_labels)),
    k_fill = dplyr::if_else(scheme == scheme_labels[["k6"]], 6L, 5L)
  )

p2 <- ggplot(prof_long, aes(x = guild, y = mean_rel, fill = amd_clusts)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = cts_cols_6, name = "amd_clusts", drop = FALSE) +
  facet_wrap(~scheme, ncol = 1) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "A7: Guild profiles by CTS (CTS1 = highest euplanctonic)",
    subtitle = paste(
      "Manuscript merge: rank k=6 raw clusters by mean total_nspp_by_lake_core;",
      "merge richness ranks 6 into 5 (main_script.Rmd); then euplanctonic CTS labels"
    ),
    y = "Mean relative abundance"
  )
ggsave("outputs/revision/figures/A7_cts_guild_profiles.png", p2, width = 10, height = 9, dpi = 150)
ggsave(
  "outputs/revision/figures/A7_cts_guild_profiles_k5_k6_merged.png",
  p2,
  width = 10,
  height = 9,
  dpi = 150
)

message("A7 complete (three schemes; manuscript k=5 = merged from shared k=6 partition)")
