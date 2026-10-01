# A6 — Volcanism sensitivity: exclude ±1 thirty-year bin around tephra ages
source("scripts/revision/_bootstrap.R")
shared <- revision_bootstrap()

message("=== A6 volcanism sensitivity ===")

tephra <- read.csv("data/revision/tephra_events.csv", stringsAsFactors = FALSE) %>%
  dplyr::mutate(age_ce = as.numeric(age_ce)) %>%
  dplyr::filter(is.finite(age_ce))

bin_size <- 30
tephra_bins <- tephra %>%
  dplyr::mutate(
    age_bin = floor(age_ce / bin_size) * bin_size,
    exclude_lo = age_bin - bin_size,
    exclude_hi = age_bin + bin_size
  )
revision_write_csv(tephra_bins, "outputs/revision/A6_tephra_bins.csv")

exclude_ages <- sort(unique(c(
  tephra_bins$age_bin,
  tephra_bins$exclude_lo,
  tephra_bins$exclude_hi
)))
exclude_ages <- exclude_ages[is.finite(exclude_ages)]

joined <- shared$joined_df_30yr %>%
  dplyr::mutate(tephra_flag = age_ce %in% exclude_ages)

flagged <- joined %>% dplyr::filter(tephra_flag)
revision_write_csv(flagged %>% dplyr::select(age_ce, n_lakes, tephra_flag),
                   "outputs/revision/A6_flagged_bins.csv")

joined_ex <- joined %>% dplyr::filter(!tephra_flag)

guild_cols <- revision_guild_cols
env_predictors <- c("NAO_Median_Value", "estimate")

vp_full <- calc_phase_varpart(
  joined_df = shared$joined_df_30yr,
  guild_cols = guild_cols,
  env_predictors = env_predictors,
  phase_breaks = revision_phase_breaks,
  min_n = 5
) %>% dplyr::mutate(scenario = "including_tephra_bins")

vp_ex <- calc_phase_varpart(
  joined_df = joined_ex,
  guild_cols = guild_cols,
  env_predictors = env_predictors,
  phase_breaks = revision_phase_breaks,
  min_n = 5
) %>% dplyr::mutate(scenario = "excluding_tephra_pm1bin")

vp <- dplyr::bind_rows(vp_full, vp_ex)
revision_write_csv(vp, "outputs/revision/A6_varpart_tephra_sensitivity.csv")

delta <- vp_full %>%
  dplyr::select(phase, component, value_full = value) %>%
  dplyr::left_join(
    vp_ex %>% dplyr::select(phase, component, value_ex = value),
    by = c("phase", "component")
  ) %>%
  dplyr::mutate(delta = value_ex - value_full)
revision_write_csv(delta, "outputs/revision/A6_varpart_delta.csv")
print(delta)

# Annotation table for Fig. 2 overlay
annot <- tephra_bins %>%
  dplyr::distinct(age_ce, event, source) %>%
  dplyr::arrange(age_ce)
revision_write_csv(annot, "outputs/revision/A6_fig2_eruption_annotations.csv")

p <- ggplot(vp %>% dplyr::filter(!is.na(value)),
            aes(x = phase, y = value, fill = component)) +
  geom_col(position = "stack") +
  facet_wrap(~scenario) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "A6: Varpart with/without tephra-neighbour bins", y = "Adj. R²")
ggsave("outputs/revision/figures/A6_tephra_varpart.png", p, width = 9, height = 5, dpi = 150)

message("A6 complete (verify tephra ages in data/revision/tephra_events.csv against Fig. S5)")
