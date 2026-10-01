# A10 — Fig. 3e analogue: observed vs rarefied species richness by CTS
# Side-by-side boxplots for each Community Trophic Structure.
source("scripts/revision/_bootstrap.R")
shared <- revision_bootstrap()

message("=== A10 Fig 3e observed vs rarefied richness by CTS ===")

suppressPackageStartupMessages({
  library(multcompView)
  library(ggpubr)
})

core_lakes <- c(
  "Azul", "Caldeirao", "Caveiro", "Empadadas Norte",
  "Funda", "Ginjal", "Peixinho", "Prata", "Santiago"
)

# ---- CTS assignments (k = 5, matching manuscript Fig. 3) ----
cts_path <- "outputs/revision/A7_sample_cts5_assignments.csv"
if (!file.exists(cts_path)) {
  message("CTS assignments missing; running A7...")
  source("scripts/revision/A07_split_cts5_cts6.R")
}
cts <- readr::read_csv(cts_path, show_col_types = FALSE) %>%
  dplyr::filter(lake %in% core_lakes, is.finite(age_ce), age_ce > 0) %>%
  dplyr::mutate(amd_clusts = as.character(amd_clusts))

# ---- Sample-level observed + rarefied richness (diatoms + chironomids) ----
diat <- revision_bind_wide_with_ages(
  shared$ls_df_diat_wide_codes[intersect(names(shared$ls_df_diat_wide_codes), core_lakes)],
  shared$ls_df_codes_diat[intersect(names(shared$ls_df_codes_diat), core_lakes)]
)
chiro <- revision_bind_wide_with_ages(
  shared$ls_df_chiro_wide_codes[intersect(names(shared$ls_df_chiro_wide_codes), core_lakes)],
  shared$df_chiro_codes[intersect(names(shared$df_chiro_codes), core_lakes)]
)

tax_d <- revision_taxon_cols(diat)
tax_c <- revision_taxon_cols(chiro)

diat_n <- revision_choose_rarefy_n(
  diat$total_count,
  target_keep = 0.70,
  candidates = c(100L, 200L, 300L, 400L, 500L)
)
chiro_n <- revision_choose_rarefy_n(
  chiro$total_count,
  target_keep = 0.70,
  candidates = c(5L, 10L, 15L, 20L, 25L, 30L, 40L, 50L)
)

rarefy_choice <- tibble::tibble(
  group = c("diatoms", "chironomids"),
  rarefy_n = c(diat_n, chiro_n),
  n_samples = c(nrow(diat), nrow(chiro)),
  n_kept = c(
    sum(diat$total_count >= diat_n, na.rm = TRUE),
    sum(chiro$total_count >= chiro_n, na.rm = TRUE)
  ),
  frac_kept = c(
    mean(diat$total_count >= diat_n, na.rm = TRUE),
    mean(chiro$total_count >= chiro_n, na.rm = TRUE)
  )
)
revision_write_csv(rarefy_choice, "outputs/revision/A10_sample_rarefy_depths.csv")
print(rarefy_choice)

message("Rarefying diatom samples to n=", diat_n, " ...")
diat_rich <- purrr::map_dfr(seq_len(nrow(diat)), function(i) {
  counts <- as.numeric(diat[i, tax_d])
  tibble::tibble(
    lake = diat$lake[i],
    core_depth_id = diat$core_depth_id[i],
    diat_observed = sum(counts > 0, na.rm = TRUE),
    diat_count = diat$total_count[i],
    diat_rarefied = revision_rarefy_richness(counts, n = diat_n, n_rep = 10L, seed = i)
  )
})

message("Rarefying chironomid samples to n=", chiro_n, " ...")
chiro_rich <- purrr::map_dfr(seq_len(nrow(chiro)), function(i) {
  counts <- as.numeric(chiro[i, tax_c])
  tibble::tibble(
    lake = chiro$lake[i],
    core_depth_id = chiro$core_depth_id[i],
    chiro_observed = sum(counts > 0, na.rm = TRUE),
    chiro_count = chiro$total_count[i],
    chiro_rarefied = revision_rarefy_richness(counts, n = chiro_n, n_rep = 10L, seed = i)
  )
})

rich <- diat_rich %>%
  dplyr::inner_join(chiro_rich, by = c("lake", "core_depth_id")) %>%
  dplyr::mutate(
    observed_richness = diat_observed + chiro_observed,
    rarefied_richness = dplyr::if_else(
      is.finite(diat_rarefied) & is.finite(chiro_rarefied),
      diat_rarefied + chiro_rarefied,
      NA_real_
    )
  )

# Prefer stored total_nspp when available (matches Fig. 3e exactly)
if (is.null(shared$df_fgroups_glob_div)) {
  e <- new.env(parent = emptyenv())
  load("data/clean_source_data_files.RData", envir = e)
  if (exists("df_fgroups_glob_div", envir = e, inherits = FALSE)) {
    shared$df_fgroups_glob_div <- get("df_fgroups_glob_div", envir = e)
  }
}
stored_totals <- shared$df_fgroups_glob_div
if (!is.null(stored_totals)) {
  stored <- stored_totals %>%
    dplyr::distinct(lake, core_depth_id, total_nspp_by_lake_core) %>%
    dplyr::filter(is.finite(total_nspp_by_lake_core))
  rich <- rich %>%
    dplyr::left_join(stored, by = c("lake", "core_depth_id")) %>%
    dplyr::mutate(
      observed_richness = dplyr::coalesce(
        as.numeric(total_nspp_by_lake_core),
        observed_richness
      )
    )
}

plot_df <- cts %>%
  dplyr::inner_join(
    rich %>% dplyr::select(lake, core_depth_id, observed_richness, rarefied_richness,
                           diat_observed, chiro_observed, diat_rarefied, chiro_rarefied,
                           diat_count, chiro_count),
    by = c("lake", "core_depth_id")
  ) %>%
  dplyr::filter(is.finite(observed_richness))

# Re-label CTS by mean observed richness (manuscript Fig. 3 convention:
# CTS1 = lowest richness … CTS5 = highest), then plot CTS5→CTS1 left→right.
rank_map <- plot_df %>%
  dplyr::group_by(amd_clusts) %>%
  dplyr::summarise(mean_obs = mean(observed_richness, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(mean_obs) %>%
  dplyr::mutate(cts_ranked = paste0("CTS", dplyr::row_number()))

revision_write_csv(rank_map, "outputs/revision/A10_cts_richness_rank_map.csv")

plot_df <- plot_df %>%
  dplyr::left_join(rank_map %>% dplyr::select(amd_clusts, cts_ranked), by = "amd_clusts") %>%
  dplyr::mutate(amd_clusts = cts_ranked) %>%
  dplyr::select(-cts_ranked)

# Order CTS5 → CTS1 as in Fig. 3e
cts_levels <- paste0("CTS", 5:1)
plot_df <- plot_df %>%
  dplyr::mutate(amd_clusts = factor(amd_clusts, levels = cts_levels))

long_df <- plot_df %>%
  tidyr::pivot_longer(
    cols = c(observed_richness, rarefied_richness),
    names_to = "richness_type",
    values_to = "n_species"
  ) %>%
  dplyr::mutate(
    richness_type = dplyr::recode(
      richness_type,
      observed_richness = "Observed",
      rarefied_richness = "Rarefied"
    ),
    richness_type = factor(richness_type, levels = c("Observed", "Rarefied"))
  ) %>%
  dplyr::filter(is.finite(n_species))

revision_write_csv(
  plot_df %>%
    dplyr::select(
      lake, core_depth_id, age_ce, amd_clusts,
      observed_richness, rarefied_richness,
      diat_observed, chiro_observed, diat_rarefied, chiro_rarefied,
      diat_count, chiro_count
    ),
  "outputs/revision/A10_cts_observed_vs_rarefied_richness.csv"
)

# Summary stats
summ <- long_df %>%
  dplyr::group_by(amd_clusts, richness_type) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean = mean(n_species),
    median = stats::median(n_species),
    sd = stats::sd(n_species),
    .groups = "drop"
  )
revision_write_csv(summ, "outputs/revision/A10_cts_richness_summary.csv")
print(summ)

# Letters for observed and rarefied separately (Tukey)
make_letters <- function(df, ycol, type_label) {
  d <- df %>%
    dplyr::filter(richness_type == type_label, is.finite(n_species)) %>%
    dplyr::mutate(amd_clusts = droplevels(amd_clusts))
  if (dplyr::n_distinct(d$amd_clusts) < 2 || nrow(d) < 5) {
    return(tibble::tibble(
      amd_clusts = levels(df$amd_clusts),
      letters = NA_character_,
      richness_type = type_label,
      y_position = NA_real_
    ))
  }
  fit <- stats::aov(n_species ~ amd_clusts, data = d)
  tuk <- stats::TukeyHSD(fit)
  cld <- multcompView::multcompLetters4(fit, tuk)
  tibble::tibble(
    amd_clusts = names(cld$amd_clusts$Letters),
    letters = unname(cld$amd_clusts$Letters),
    richness_type = type_label
  ) %>%
    dplyr::left_join(
      d %>% dplyr::group_by(amd_clusts) %>% dplyr::summarise(max_y = max(n_species), .groups = "drop"),
      by = "amd_clusts"
    ) %>%
    dplyr::mutate(
      amd_clusts = factor(amd_clusts, levels = cts_levels),
      richness_type = factor(richness_type, levels = c("Observed", "Rarefied")),
      y_position = max_y * 1.08
    )
}

letters_df <- dplyr::bind_rows(
  make_letters(long_df, "n_species", "Observed"),
  make_letters(long_df, "n_species", "Rarefied")
)

anova_obs <- summary(stats::aov(n_species ~ amd_clusts, data = dplyr::filter(long_df, richness_type == "Observed")))
anova_rar <- summary(stats::aov(n_species ~ amd_clusts, data = dplyr::filter(long_df, richness_type == "Rarefied")))
p_obs <- anova_obs[[1]][["Pr(>F)"]][1]
p_rar <- anova_rar[[1]][["Pr(>F)"]][1]

anova_tbl <- tibble::tibble(
  richness_type = c("Observed", "Rarefied"),
  anova_p = c(p_obs, p_rar),
  diat_rarefy_n = diat_n,
  chiro_rarefy_n = chiro_n
)
revision_write_csv(anova_tbl, "outputs/revision/A10_anova_pvalues.csv")

# Primary figure: CTS colours for Observed; lighter companion for Rarefied
obs_cols <- setNames(viridis::viridis(5, direction = 1), cts_levels)
# Build paired palette: Observed = full CTS colour; Rarefied = desaturated companion
pair_levels <- interaction(
  factor(rep(cts_levels, each = 2), levels = cts_levels),
  factor(rep(c("Observed", "Rarefied"), times = 5), levels = c("Observed", "Rarefied")),
  sep = " · ",
  lex.order = TRUE
)

long_df <- long_df %>%
  dplyr::mutate(
    pair = interaction(amd_clusts, richness_type, sep = " · ", lex.order = TRUE)
  )

pair_cols <- unlist(lapply(cts_levels, function(cts) {
  base <- grDevices::col2rgb(obs_cols[[cts]]) / 255
  obs <- obs_cols[[cts]]
  rar <- grDevices::rgb(
    0.55 * base[1] + 0.45,
    0.55 * base[2] + 0.45,
    0.55 * base[3] + 0.45
  )
  stats::setNames(c(obs, rar), paste(cts, c("Observed", "Rarefied"), sep = " · "))
}))

p <- ggplot(long_df, aes(x = amd_clusts, y = n_species, fill = pair)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.7,
    outlier.size = 0.7,
    colour = "grey20",
    linewidth = 0.3
  ) +
  geom_text(
    data = letters_df %>% dplyr::filter(is.finite(y_position)),
    aes(x = amd_clusts, y = y_position, label = letters, group = richness_type),
    position = position_dodge(width = 0.8),
    inherit.aes = FALSE,
    size = 3.5
  ) +
  scale_fill_manual(values = pair_cols, guide = "none") +
  scale_x_discrete(limits = cts_levels) +
  annotate(
    "text",
    x = 0.6,
    y = max(long_df$n_species, na.rm = TRUE) * 0.08,
    label = paste0(
      "Observed ANOVA p ", ifelse(p_obs < 2.2e-16, "< 2.2e-16", paste0("= ", signif(p_obs, 3))),
      "\nRarefied ANOVA p ", ifelse(p_rar < 2.2e-16, "< 2.2e-16", paste0("= ", signif(p_rar, 3)))
    ),
    hjust = 0,
    size = 3.2
  ) +
  labs(
    title = "Fig. 3e revision: observed vs rarefied species richness by CTS",
    subtitle = paste0(
      "Left bar = observed; right bar = rarefied ",
      "(diatoms n=", diat_n, "; chironomids n=", chiro_n, ")"
    ),
    x = "Community Trophic Structures",
    y = "Number of species"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(colour = "grey70", linewidth = 0.2)
  )

ggsave(
  "outputs/revision/figures/A10_fig3e_observed_vs_rarefied_cts.png",
  p,
  width = 8.5,
  height = 5.5,
  dpi = 200
)

# Legend-friendly companion (Observed vs Rarefied only)
p2 <- ggplot(long_df, aes(x = amd_clusts, y = n_species, fill = richness_type)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.65, outlier.size = 0.7) +
  geom_text(
    data = letters_df %>% dplyr::filter(is.finite(y_position)),
    aes(x = amd_clusts, y = y_position, label = letters, group = richness_type),
    position = position_dodge(width = 0.75),
    inherit.aes = FALSE,
    size = 3.5
  ) +
  scale_fill_manual(
    name = NULL,
    values = c(Observed = "#440154", Rarefied = "#35B779"),
    labels = c(
      Observed = "Observed",
      Rarefied = paste0("Rarefied (diat n=", diat_n, "; chiro n=", chiro_n, ")")
    )
  ) +
  scale_x_discrete(limits = cts_levels) +
  labs(
    title = "Fig. 3e revision: observed vs rarefied richness by CTS",
    subtitle = "Paired boxplots per CTS (Tukey letters above each series)",
    x = "Community Trophic Structures",
    y = "Number of species"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  "outputs/revision/figures/A10_fig3e_observed_vs_rarefied_cts_dodged.png",
  p2,
  width = 8.5,
  height = 5.5,
  dpi = 200
)

message("A10 complete")
