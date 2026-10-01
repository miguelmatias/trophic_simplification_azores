# A7 — Split CTS5/CTS6 (compare k=5 lumped vs k=6 unlumped)
source("scripts/revision/_bootstrap.R")
shared <- revision_bootstrap()

message("=== A7 CTS5/CTS6 split ===")

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
  replace(is.na(.), 0)

set.seed(42)
message("Clustering k=6 (no lumping)...")
cts6 <- revision_assign_cts(norm_abund, k = 6L, iterations = 300L)
message("Clustering k=5...")
cts5 <- revision_assign_cts(norm_abund, k = 5L, iterations = 300L)

# Also create lumped version of k=6 (CTS6 -> CTS5) to mimic manuscript decision
cts6_lumped <- cts6 %>%
  dplyr::mutate(
    amd_clusts = as.character(amd_clusts),
    amd_clusts = ifelse(amd_clusts == "CTS6", "CTS5", amd_clusts),
    amd_clusts = factor(amd_clusts, levels = paste0("CTS", 1:5))
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
  profiles(cts5, "k5"),
  profiles(cts6, "k6"),
  profiles(cts6_lumped, "k6_lumped_to_5")
)
revision_write_csv(prof, "outputs/revision/A7_cts_guild_profiles.csv")

occ <- dplyr::bind_rows(
  revision_cts_occupancy(cts5) %>% dplyr::mutate(scheme = "k5"),
  revision_cts_occupancy(cts6) %>% dplyr::mutate(scheme = "k6")
)
revision_write_csv(occ, "outputs/revision/A7_cts_occupancy_by_lake.csv")

# Temporal occupancy (regional proportion of samples per CTS in 30-yr bins)
temporal <- function(df, label) {
  df %>%
    dplyr::filter(age_ce > 0) %>%
    dplyr::mutate(age_bin = floor(age_ce / 30) * 30) %>%
    dplyr::count(age_bin, amd_clusts, name = "n") %>%
    dplyr::group_by(age_bin) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(scheme = label) %>%
    dplyr::rename(age_ce = age_bin)
}

temp <- dplyr::bind_rows(temporal(cts5, "k5"), temporal(cts6, "k6"))
revision_write_csv(temp, "outputs/revision/A7_cts_temporal_proportions.csv")

# Early-phase stability: fraction of samples in "simple" CTS classes before 750 CE
early_stability <- dplyr::bind_rows(
  cts5 %>% dplyr::mutate(scheme = "k5"),
  cts6 %>% dplyr::mutate(scheme = "k6")
) %>%
  dplyr::filter(age_ce > 0, age_ce < 750) %>%
  dplyr::count(scheme, amd_clusts) %>%
  dplyr::group_by(scheme) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup()
revision_write_csv(early_stability, "outputs/revision/A7_early_phase_cts_composition.csv")
print(early_stability)

# Save sample assignments for downstream A8
revision_write_csv(
  cts6 %>% dplyr::select(lake, core_depth_id, age_ce, amd_clusts),
  "outputs/revision/A7_sample_cts6_assignments.csv"
)
revision_write_csv(
  cts5 %>% dplyr::select(lake, core_depth_id, age_ce, amd_clusts),
  "outputs/revision/A7_sample_cts5_assignments.csv"
)

p <- ggplot(temp, aes(x = age_ce, y = prop, fill = amd_clusts)) +
  geom_area(alpha = 0.9) +
  facet_wrap(~scheme, ncol = 1) +
  theme_minimal(base_size = 11) +
  labs(title = "A7: CTS temporal occupancy (k=5 vs k=6)", x = "Age (CE)", y = "Proportion")
ggsave("outputs/revision/figures/A7_cts_temporal.png", p, width = 9, height = 7, dpi = 150)

p2 <- prof %>%
  dplyr::filter(scheme %in% c("k5", "k6")) %>%
  tidyr::pivot_longer(dplyr::all_of(revision_guild_cols), names_to = "guild", values_to = "mean_rel") %>%
  ggplot(aes(x = guild, y = mean_rel, fill = amd_clusts)) +
  geom_col(position = "dodge") +
  facet_wrap(~scheme, ncol = 1) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "A7: Guild profiles by CTS", y = "Mean relative abundance")
ggsave("outputs/revision/figures/A7_cts_guild_profiles.png", p2, width = 10, height = 7, dpi = 150)

message("A7 complete")
