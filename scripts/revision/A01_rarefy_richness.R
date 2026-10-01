# A1 — Rarefy richness in 30-yr bins; redraw revision Fig 3e / Fig 4 style plots
source("scripts/revision/_bootstrap.R")
shared <- revision_bootstrap()

message("=== A1 rarefaction ===")

core_lakes <- c(
  "Azul", "Caldeirao", "Caveiro", "Empadadas Norte",
  "Funda", "Ginjal", "Peixinho", "Prata", "Santiago"
)

# ---- Sample-level count diagnostics ----
diat_ages <- revision_bind_wide_with_ages(
  shared$ls_df_diat_wide_codes[intersect(names(shared$ls_df_diat_wide_codes), core_lakes)],
  shared$ls_df_codes_diat[intersect(names(shared$ls_df_codes_diat), core_lakes)]
)
chiro_ages <- revision_bind_wide_with_ages(
  shared$ls_df_chiro_wide_codes[intersect(names(shared$ls_df_chiro_wide_codes), core_lakes)],
  shared$df_chiro_codes[intersect(names(shared$df_chiro_codes), core_lakes)]
)

sample_counts <- dplyr::bind_rows(
  diat_ages %>%
    dplyr::select(lake, core_depth_id, age_ce, total_count) %>%
    dplyr::mutate(group = "Producers"),
  chiro_ages %>%
    dplyr::select(lake, core_depth_id, age_ce, total_count) %>%
    dplyr::mutate(group = "Consumers")
)
revision_write_csv(sample_counts, "outputs/revision/A1_sample_counts.csv")

count_summary <- sample_counts %>%
  dplyr::group_by(group, lake) %>%
  dplyr::summarise(
    n_samples = dplyr::n(),
    min_count = min(total_count, na.rm = TRUE),
    median_count = stats::median(total_count, na.rm = TRUE),
    mean_count = mean(total_count, na.rm = TRUE),
    max_count = max(total_count, na.rm = TRUE),
    .groups = "drop"
  )
revision_write_csv(count_summary, "outputs/revision/A1_count_summary_by_lake.csv")

# ---- Pool into 30-yr bins ----
diat_bins <- revision_pool_bins(diat_ages)
chiro_bins <- revision_pool_bins(chiro_ages)

n_prod <- revision_choose_rarefy_n(diat_bins$total_count, target_keep = 0.70)
n_cons <- revision_choose_rarefy_n(
  chiro_bins$total_count,
  target_keep = 0.70,
  candidates = c(5L, 10L, 15L, 20L, 25L, 30L, 40L, 50L, 75L, 100L)
)

rarefy_choice <- tibble::tibble(
  group = c("Producers", "Consumers"),
  rarefy_n = c(n_prod, n_cons),
  n_bins_total = c(nrow(diat_bins), nrow(chiro_bins)),
  n_bins_kept = c(
    sum(diat_bins$total_count >= n_prod, na.rm = TRUE),
    sum(chiro_bins$total_count >= n_cons, na.rm = TRUE)
  ),
  frac_bins_kept = c(
    mean(diat_bins$total_count >= n_prod, na.rm = TRUE),
    mean(chiro_bins$total_count >= n_cons, na.rm = TRUE)
  )
)
revision_write_csv(rarefy_choice, "outputs/revision/A1_rarefy_depth_choice.csv")
print(rarefy_choice)

rarefy_one_group <- function(bins, n, group_label) {
  tax <- revision_taxon_cols(bins)
  purrr::pmap_dfr(
    list(seq_len(nrow(bins))),
    function(i) {
      row <- bins[i, ]
      richness <- revision_rarefy_richness(as.numeric(row[tax]), n = n, n_rep = 10L, seed = i)
      tibble::tibble(
        lake = row$lake,
        age_ce = row$age_ce,
        group = group_label,
        n_samples_in_bin = row$n_samples,
        total_count = row$total_count,
        rarefy_n = n,
        rarefied_richness = richness,
        observed_richness = sum(as.numeric(row[tax]) > 0)
      )
    }
  )
}

message("Rarefying producers to n=", n_prod, " ...")
prod_rare <- rarefy_one_group(diat_bins, n_prod, "Producers")
message("Rarefying consumers to n=", n_cons, " ...")
cons_rare <- rarefy_one_group(chiro_bins, n_cons, "Consumers")

rare_long <- dplyr::bind_rows(prod_rare, cons_rare)
revision_write_csv(rare_long, "outputs/revision/A1_rarefied_richness_bins.csv")

# Regional mean (no within-lake min-max)
regional <- rare_long %>%
  dplyr::filter(is.finite(rarefied_richness)) %>%
  dplyr::group_by(group, age_ce) %>%
  dplyr::summarise(
    n_lakes = dplyr::n_distinct(lake),
    mean_richness = mean(rarefied_richness, na.rm = TRUE),
    se_richness = stats::sd(rarefied_richness, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  )
revision_write_csv(regional, "outputs/revision/A1_regional_rarefied_richness.csv")

# Sensitivity: alternate rarefaction depths (lighter reps)
sens_depths_prod <- unique(c(300L, 400L, n_prod))
sens_depths_cons <- unique(c(20L, 30L, n_cons))
rarefy_one_group_light <- function(bins, n, group_label) {
  tax <- revision_taxon_cols(bins)
  purrr::map_dfr(seq_len(nrow(bins)), function(i) {
    row <- bins[i, ]
    richness <- revision_rarefy_richness(as.numeric(row[tax]), n = n, n_rep = 5L, seed = i)
    tibble::tibble(
      lake = row$lake, age_ce = row$age_ce, group = group_label,
      total_count = row$total_count, rarefy_n = n, rarefied_richness = richness
    )
  })
}
sens_prod <- purrr::map_dfr(sens_depths_prod, function(n) {
  rarefy_one_group_light(diat_bins, n, "Producers") %>% dplyr::mutate(sensitivity_n = n)
})
sens_cons <- purrr::map_dfr(sens_depths_cons, function(n) {
  rarefy_one_group_light(chiro_bins, n, "Consumers") %>% dplyr::mutate(sensitivity_n = n)
})
sens <- dplyr::bind_rows(sens_prod, sens_cons)
revision_write_csv(sens, "outputs/revision/A1_rarefaction_sensitivity.csv")

# Trend summaries (pre/post 1600)
trend_summary <- rare_long %>%
  dplyr::filter(is.finite(rarefied_richness)) %>%
  dplyr::mutate(period = ifelse(age_ce < 1600, "pre_1600", "post_1600")) %>%
  dplyr::group_by(group, period) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_richness = mean(rarefied_richness),
    median_richness = stats::median(rarefied_richness),
    .groups = "drop"
  )
revision_write_csv(trend_summary, "outputs/revision/A1_period_richness_summary.csv")

# Figures (revision versions)
p_counts <- ggplot(sample_counts, aes(x = lake, y = total_count, fill = group)) +
  geom_boxplot(outlier.size = 0.6, alpha = 0.8) +
  facet_wrap(~group, scales = "free_y") +
  coord_flip() +
  theme_minimal(base_size = 11) +
  labs(title = "A1: Count totals per sample", y = "Count", x = NULL) +
  theme(legend.position = "none")
ggsave("outputs/revision/figures/A1_sample_counts.png", p_counts, width = 9, height = 6, dpi = 150)

p_fig4 <- ggplot(regional, aes(x = age_ce, y = mean_richness, colour = group, fill = group)) +
  geom_ribbon(aes(ymin = mean_richness - se_richness, ymax = mean_richness + se_richness), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~group, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 11) +
  labs(
    title = "A1: Regional rarefied richness (30-yr bins; no min–max scaling)",
    subtitle = paste0("Producers n=", n_prod, "; Consumers n=", n_cons),
    x = "Age (CE)", y = "Rarefied richness"
  ) +
  theme(legend.position = "none")
ggsave("outputs/revision/figures/A1_fig4_rarefied_richness.png", p_fig4, width = 8, height = 7, dpi = 150)

# Fig 3e style: lake trajectories of rarefied richness
p_3e <- ggplot(rare_long %>% dplyr::filter(is.finite(rarefied_richness)),
               aes(x = age_ce, y = rarefied_richness, colour = group)) +
  geom_line(alpha = 0.85) +
  facet_wrap(~lake, scales = "free_y") +
  theme_minimal(base_size = 10) +
  labs(
    title = "A1: Lake-level rarefied richness (revision Fig. 3e analogue)",
    x = "Age (CE)", y = "Rarefied richness"
  )
ggsave("outputs/revision/figures/A1_fig3e_lake_rarefied_richness.png", p_3e, width = 11, height = 8, dpi = 150)

message("A1 complete")
