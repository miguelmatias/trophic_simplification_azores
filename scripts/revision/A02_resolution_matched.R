# A2 — Resolution-matched producer/consumer comparison
source("scripts/revision/_bootstrap.R")
shared <- revision_bootstrap()

message("=== A2 resolution matching ===")

core_lakes <- c(
  "Azul", "Caldeirao", "Caveiro", "Empadadas Norte",
  "Funda", "Ginjal", "Peixinho", "Prata", "Santiago"
)

diat_wide <- shared$ls_df_diat_wide_codes[core_lakes]
diat_codes <- shared$ls_df_codes_diat[core_lakes]
chiro_wide <- shared$ls_df_chiro_wide_codes[core_lakes]
chiro_codes <- shared$df_chiro_codes[core_lakes]

# Species-level DCA (match manuscript convention: flip Prata diatoms)
message("Species-level DCA...")
diat_sp <- revision_compute_dca1_list(diat_wide, diat_codes, flip_lakes = "Prata")
chiro_sp <- revision_compute_dca1_list(chiro_wide, chiro_codes)

# Genus-level diatoms
message("Genus aggregation + DCA...")
diat_genus_wide <- revision_aggregate_diatoms_to_genus(diat_wide)
# ensure lake column present
diat_genus_wide <- purrr::imap(diat_genus_wide, function(df, nm) {
  if (!"lake" %in% names(df) || all(is.na(df$lake))) df$lake <- nm
  df
})
diat_gen <- revision_compute_dca1_list(diat_genus_wide, diat_codes, flip_lakes = "Prata")

# Guild-space DCA (4 producer + 4 consumer guilds separately)
message("Guild-space DCA...")
guild_prod_cols <- c("high_profile", "low_profile", "motile", "euplanktonic")
guild_cons_cols <- c("algivore", "detritivore", "plantivore", "predator")

fg <- shared$df_fgroups %>%
  dplyr::filter(lake %in% core_lakes, is.finite(age_ce)) %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(c(guild_prod_cols, guild_cons_cols)), ~ replace(.x, is.na(.x), 0)))

guild_dca_one <- function(df, cols, group_label) {
  purrr::map_dfr(split(df, df$lake), function(lk_df) {
    mat <- lk_df %>%
      dplyr::select(dplyr::all_of(cols)) %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), ~ .x / pmax(rowSums(dplyr::across(dplyr::everything())), 1e-8)))
    rs <- rowSums(mat)
    keep_rows <- rs > 0 & is.finite(rs)
    if (sum(keep_rows) < 3) {
      return(tibble::tibble(
        lake = lk_df$lake, core_depth_id = lk_df$core_depth_id,
        age_ce = lk_df$age_ce, DCA1 = NA_real_, group = group_label
      ))
    }
    mat <- mat[keep_rows, , drop = FALSE]
    meta <- lk_df[keep_rows, , drop = FALSE]
    hell <- vegan::decostand(mat, method = "hellinger")
    keep <- vapply(hell, function(x) stats::sd(x) > 0, logical(1))
    if (sum(keep) < 2) {
      return(tibble::tibble(
        lake = meta$lake, core_depth_id = meta$core_depth_id,
        age_ce = meta$age_ce, DCA1 = NA_real_, group = group_label
      ))
    }
    dca <- vegan::decorana(hell[, keep, drop = FALSE])
    sc <- as.data.frame(vegan::scores(dca, display = "sites", choices = 1))
    names(sc) <- "DCA1"
    tibble::tibble(
      lake = meta$lake,
      core_depth_id = meta$core_depth_id,
      age_ce = meta$age_ce,
      DCA1 = as.numeric(sc$DCA1),
      group = group_label
    )
  })
}

guild_prod <- guild_dca_one(fg, guild_prod_cols, "Producers")
guild_cons <- guild_dca_one(fg, guild_cons_cols, "Consumers")

run_timing <- function(diat_df, chiro_df, label) {
  message("Timing: ", label)
  out <- compute_lake_level_gam_timing(
    diat_df = diat_df %>% dplyr::select(lake, age_ce, DCA1),
    chiro_df = chiro_df %>% dplyr::select(lake, age_ce, DCA1),
    lakes = core_lakes,
    exclude_lakes = c("Fogo", "Furnas"),
    min_n = 8L,
    onset_min_ce = 800
  )
  out$comparison <- out$comparison %>% dplyr::mutate(resolution = label)
  out$summary <- out$summary %>% dplyr::mutate(resolution = label)
  out
}

t_sp <- run_timing(diat_sp, chiro_sp, "species")
t_gen <- run_timing(diat_gen, chiro_sp, "genus_diatoms")
t_guild <- run_timing(
  guild_prod %>% dplyr::select(lake, age_ce, DCA1),
  guild_cons %>% dplyr::select(lake, age_ce, DCA1),
  "guild"
)

comparison <- dplyr::bind_rows(t_sp$comparison, t_gen$comparison, t_guild$comparison)
summary_tbl <- dplyr::bind_rows(t_sp$summary, t_gen$summary, t_guild$summary)

revision_write_csv(comparison, "outputs/revision/A2_onset_comparison_by_resolution.csv")
revision_write_csv(summary_tbl, "outputs/revision/A2_timing_summary_by_resolution.csv")

leader_counts <- comparison %>%
  dplyr::count(resolution, onset_leader) %>%
  tidyr::pivot_wider(names_from = onset_leader, values_from = n, values_fill = 0)
revision_write_csv(leader_counts, "outputs/revision/A2_onset_leader_counts.csv")

mag_summary <- comparison %>%
  dplyr::group_by(resolution) %>%
  dplyr::summarise(
    n_lakes = dplyr::n(),
    mean_mag_diff = mean(magnitude_diff, na.rm = TRUE),
    median_mag_diff = stats::median(magnitude_diff, na.rm = TRUE),
    n_producers_stronger = sum(magnitude_diff > 0, na.rm = TRUE),
    n_consumers_stronger = sum(magnitude_diff < 0, na.rm = TRUE),
    n_producers_earlier = sum(onset_leader == "Producers earlier", na.rm = TRUE),
    n_consumers_earlier = sum(onset_leader == "Consumers earlier", na.rm = TRUE),
    .groups = "drop"
  )
revision_write_csv(mag_summary, "outputs/revision/A2_magnitude_onset_summary.csv")
print(mag_summary)

p <- ggplot(comparison, aes(x = resolution, y = magnitude_diff, colour = onset_leader)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_jitter(width = 0.1, height = 0, size = 2.5) +
  theme_minimal(base_size = 11) +
  labs(
    title = "A2: Producer−consumer DCA change difference by resolution",
    y = "total_dca_change_P − total_dca_change_C",
    x = NULL
  )
ggsave("outputs/revision/figures/A2_magnitude_by_resolution.png", p, width = 8, height = 5, dpi = 150)

message("A2 complete")
