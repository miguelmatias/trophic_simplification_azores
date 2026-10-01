# A5 — Separate human vs climatic vegetation change
source("scripts/revision/_bootstrap.R")
shared <- revision_bootstrap()

message("=== A5 human vs climatic vegetation ===")

# Regional arboreal mean on 30-yr bins
arb_30 <- shared$arboreal %>%
  dplyr::mutate(age_bin = floor(age_ce / 30) * 30) %>%
  dplyr::group_by(age_bin) %>%
  dplyr::summarise(
    arboreal_pct = mean(arboreal_pct, na.rm = TRUE),
    n_sites = dplyr::n_distinct(sitename),
    .groups = "drop"
  ) %>%
  dplyr::rename(age_ce = age_bin)

arb_nao <- arb_30 %>%
  dplyr::left_join(shared$nao_30, by = "age_ce") %>%
  dplyr::left_join(shared$veg_30 %>% dplyr::select(age_ce, estimate), by = "age_ce") %>%
  dplyr::mutate(
    era = ifelse(age_ce < 1300, "pre_1300", "post_1300"),
    arboreal_loss = -arboreal_pct
  )

# Correlate arboreal % and vegetation DCA with NAO pre/post 1300
cors <- purrr::map_dfr(c("arboreal_pct", "estimate"), function(resp) {
  purrr::map_dfr(c("pre_1300", "post_1300", "all"), function(era) {
    sub <- if (era == "all") arb_nao else dplyr::filter(arb_nao, era == .env$era)
    dplyr::bind_rows(
      revision_cor_pair(sub$NAO_Median_Value, sub[[resp]], "pearson"),
      revision_cor_pair(sub$NAO_Median_Value, sub[[resp]], "spearman")
    ) %>%
      dplyr::mutate(response = resp, era = era)
  })
})
revision_write_csv(cors, "outputs/revision/A5_arboreal_nao_correlations.csv")

# Indicator onsets
ind <- shared$indicators
revision_write_csv(ind, "outputs/revision/A5_indicator_onsets_by_site.csv")

ind_summary <- ind %>%
  dplyr::group_by(category) %>%
  dplyr::summarise(
    n_sites = dplyr::n(),
    median_first_ce = stats::median(first_presence_ce, na.rm = TRUE),
    min_first_ce = min(first_presence_ce, na.rm = TRUE),
    max_first_ce = max(first_presence_ce, na.rm = TRUE),
    .groups = "drop"
  )
revision_write_csv(ind_summary, "outputs/revision/A5_indicator_onset_summary.csv")
print(ind_summary)

# Mark whether vegetation smooth change accelerates after indicator onsets
veg <- shared$veg_yearly %>%
  dplyr::mutate(
    d_est = c(NA, diff(estimate)),
    after_1300 = age_ce >= 1300
  )
veg_change <- veg %>%
  dplyr::filter(is.finite(d_est)) %>%
  dplyr::group_by(after_1300) %>%
  dplyr::summarise(
    mean_abs_change = mean(abs(d_est), na.rm = TRUE),
    mean_change = mean(d_est, na.rm = TRUE),
    .groups = "drop"
  )
revision_write_csv(veg_change, "outputs/revision/A5_veg_smooth_change_pre_post_1300.csv")

p <- ggplot(arb_nao, aes(x = age_ce)) +
  geom_line(aes(y = arboreal_pct), colour = "#1b9e77", linewidth = 0.9) +
  geom_vline(
    data = ind_summary,
    aes(xintercept = median_first_ce, colour = category),
    linetype = 2
  ) +
  theme_minimal(base_size = 11) +
  labs(
    title = "A5: Regional arboreal % with median Plantago/Cereal onsets",
    x = "Age (CE)", y = "Arboreal pollen (%)"
  )
ggsave("outputs/revision/figures/A5_arboreal_and_indicators.png", p, width = 8, height = 4.5, dpi = 150)

message("A5 complete")
