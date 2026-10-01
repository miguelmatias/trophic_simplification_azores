# A4 — NAO–vegetation correlation on 30-yr bins
source("scripts/revision/_bootstrap.R")
shared <- revision_bootstrap()

message("=== A4 NAO–vegetation correlation ===")

df <- shared$joined_df_30yr %>%
  dplyr::filter(is.finite(NAO_Median_Value), is.finite(estimate)) %>%
  dplyr::mutate(
    phase = dplyr::case_when(
      age_ce < 750 ~ "Phase 1 (<750)",
      age_ce < 1050 ~ "Phase 2 (750–1050)",
      age_ce < 1450 ~ "Phase 3 (1050–1450)",
      age_ce < 1750 ~ "Phase 4 (1450–1750)",
      TRUE ~ "Phase 5 (>1750)"
    ),
    era = ifelse(age_ce < 1300, "pre_1300", "post_1300")
  )

overall <- dplyr::bind_rows(
  revision_cor_pair(df$NAO_Median_Value, df$estimate, "pearson") %>% dplyr::mutate(slice = "overall"),
  revision_cor_pair(df$NAO_Median_Value, df$estimate, "spearman") %>% dplyr::mutate(slice = "overall")
)

by_phase <- purrr::map_dfr(split(df, df$phase), function(sub) {
  dplyr::bind_rows(
    revision_cor_pair(sub$NAO_Median_Value, sub$estimate, "pearson"),
    revision_cor_pair(sub$NAO_Median_Value, sub$estimate, "spearman")
  ) %>%
    dplyr::mutate(slice = unique(sub$phase))
})

by_era <- purrr::map_dfr(split(df, df$era), function(sub) {
  dplyr::bind_rows(
    revision_cor_pair(sub$NAO_Median_Value, sub$estimate, "pearson"),
    revision_cor_pair(sub$NAO_Median_Value, sub$estimate, "spearman")
  ) %>%
    dplyr::mutate(slice = unique(sub$era))
})

cors <- dplyr::bind_rows(overall, by_phase, by_era) %>%
  dplyr::select(slice, method, n, estimate, p.value)
revision_write_csv(cors, "outputs/revision/A4_nao_vegetation_correlations.csv")
print(cors)

p <- ggplot(df, aes(x = NAO_Median_Value, y = estimate, colour = era)) +
  geom_point(alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7) +
  facet_wrap(~era) +
  theme_minimal(base_size = 11) +
  labs(title = "A4: NAO vs vegetation change (30-yr bins)", x = "NAO", y = "Vegetation DCA1 smooth")
ggsave("outputs/revision/figures/A4_nao_vs_vegetation.png", p, width = 8, height = 4.5, dpi = 150)

message("A4 complete")
