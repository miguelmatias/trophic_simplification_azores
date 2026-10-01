# A3 — Common-period and balanced-coverage variance partitioning
source("scripts/revision/_bootstrap.R")
shared <- revision_bootstrap()

message("=== A3 coverage / common-period ===")

guild_cols <- revision_guild_cols
env_predictors <- c("NAO_Median_Value", "estimate")
phase_breaks <- revision_phase_breaks

run_varpart <- function(joined, label) {
  calc_phase_varpart(
    joined_df = joined,
    guild_cols = guild_cols,
    env_predictors = env_predictors,
    phase_breaks = phase_breaks,
    min_n = 5
  ) %>%
    dplyr::mutate(scenario = label)
}

# Coverage diagnostics
coverage <- shared$bio_lake_30 %>%
  dplyr::count(age_ce, name = "n_lakes") %>%
  dplyr::mutate(
    period = dplyr::case_when(
      age_ce < 700 ~ "pre_700",
      age_ce < 1000 ~ "700_999",
      age_ce < 1400 ~ "1000_1399",
      TRUE ~ "1400_present"
    )
  )
revision_write_csv(coverage, "outputs/revision/A3_lake_coverage_by_bin.csv")

coverage_by_period <- coverage %>%
  dplyr::group_by(period) %>%
  dplyr::summarise(
    n_bins = dplyr::n(),
    mean_n_lakes = mean(n_lakes),
    min_n_lakes = min(n_lakes),
    max_n_lakes = max(n_lakes),
    .groups = "drop"
  )
revision_write_csv(coverage_by_period, "outputs/revision/A3_coverage_by_period.csv")

# (i) full record (baseline)
vp_full <- run_varpart(shared$joined_df_30yr, "full_record")

# (i) common period 1400–2000
joined_common <- shared$joined_df_30yr %>% dplyr::filter(age_ce >= 1400, age_ce <= 2000)
vp_common <- run_varpart(joined_common, "common_1400_2000")

# (ii) bins with n_lakes >= 6
joined_n6 <- shared$joined_df_30yr %>% dplyr::filter(n_lakes >= 6)
vp_n6 <- run_varpart(joined_n6, "n_lakes_ge_6")

# (iii) island-balanced regional mean
bio_island <- revision_island_balanced_bio_reg(shared$bio_lake_30, shared$lake_island %>% dplyr::select(lake, island))
joined_island <- bio_island %>%
  dplyr::left_join(shared$nao_30, by = "age_ce") %>%
  dplyr::left_join(shared$veg_30 %>% dplyr::select(age_ce, estimate), by = "age_ce")
vp_island <- run_varpart(joined_island, "island_balanced")

# (iii-b) lake as conditioning term on sample-level 30-yr lake bins
message("Partial RDA with lake conditioning...")
lake_env <- shared$bio_lake_30 %>%
  dplyr::left_join(shared$nao_30, by = "age_ce") %>%
  dplyr::left_join(shared$veg_30 %>% dplyr::select(age_ce, estimate), by = "age_ce") %>%
  dplyr::filter(dplyr::if_all(dplyr::all_of(c(guild_cols, env_predictors)), ~ !is.na(.x)))

lake_conditioned <- purrr::imap_dfr(phase_breaks, function(bounds, phase_name) {
  lo <- bounds[1]; hi <- bounds[2]
  sub <- lake_env %>% dplyr::filter(age_ce >= lo, age_ce < hi)
  if (nrow(sub) < 8) {
    return(tibble::tibble(
      phase = phase_name, component = c("Pure NAO_Median_Value", "Pure estimate", "Shared"),
      value = NA_real_, n_bins = nrow(sub), scenario = "lake_conditioned"
    ))
  }
  comm <- as.matrix(sub[, guild_cols])
  envA <- data.frame(NAO_Median_Value = sub$NAO_Median_Value)
  envB <- data.frame(estimate = sub$estimate)
  lake_df <- data.frame(lake = factor(sub$lake))

  pure_nao <- max(0, vegan::RsquareAdj(vegan::rda(comm, envA, cbind(envB, lake_df)))$adj.r.squared)
  pure_veg <- max(0, vegan::RsquareAdj(vegan::rda(comm, envB, cbind(envA, lake_df)))$adj.r.squared)
  total <- max(0, vegan::RsquareAdj(vegan::rda(comm, cbind(envA, envB), lake_df))$adj.r.squared)
  tibble::tibble(
    phase = phase_name,
    component = c("Pure NAO_Median_Value", "Pure estimate", "Shared"),
    value = c(pure_nao, pure_veg, max(0, total - pure_nao - pure_veg)),
    n_bins = nrow(sub),
    scenario = "lake_conditioned"
  )
})

vp_all <- dplyr::bind_rows(vp_full, vp_common, vp_n6, vp_island, lake_conditioned)
revision_write_csv(vp_all, "outputs/revision/A3_varpart_by_scenario.csv")

# Wide comparison table for response letter
vp_wide <- vp_all %>%
  dplyr::mutate(component = dplyr::recode(
    component,
    "Pure NAO_Median_Value" = "Pure_NAO",
    "Pure estimate" = "Pure_vegetation"
  )) %>%
  dplyr::select(scenario, phase, component, value) %>%
  tidyr::pivot_wider(names_from = component, values_from = value)
revision_write_csv(vp_wide, "outputs/revision/A3_varpart_wide.csv")
print(vp_wide)

p <- ggplot(vp_all %>% dplyr::filter(!is.na(value)),
            aes(x = phase, y = value, fill = component)) +
  geom_col(position = "stack") +
  facet_wrap(~scenario, ncol = 2) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "A3: Phase variance partitioning under coverage scenarios", y = "Adj. R²", x = NULL)
ggsave("outputs/revision/figures/A3_varpart_scenarios.png", p, width = 10, height = 8, dpi = 150)

message("A3 complete")
