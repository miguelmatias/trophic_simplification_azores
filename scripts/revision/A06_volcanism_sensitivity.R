# A6 — Volcanism sensitivity: tephra as a predictor (not a data filter)
# Primary test: compare phase varpart community ~ NAO + vegetation
# vs community ~ NAO + vegetation + tephra (30-yr bin indicator).
# Optional diagnostic: exclusion of ±1 bin around tephra ages (not main result).
source("scripts/revision/_bootstrap.R")
shared <- revision_bootstrap()

message("=== A6 volcanism sensitivity (tephra as predictor) ===")

read_tephra <- function(path = "data/revision/tephra_events.csv") {
  if (requireNamespace("readr", quietly = TRUE)) {
    return(as.data.frame(
      readr::read_csv(path, show_col_types = FALSE,
                      locale = readr::locale(encoding = "UTF-8")),
      stringsAsFactors = FALSE
    ))
  }
  raw <- readBin(path, what = "raw", n = file.info(path)$size)
  if (length(raw) >= 3 && raw[1] == as.raw(0xef) && raw[2] == as.raw(0xbb) && raw[3] == as.raw(0xbf)) {
    raw <- raw[-(1:3)]
  }
  txt <- rawToChar(raw)
  Encoding(txt) <- "UTF-8"
  utils::read.csv(text = txt, stringsAsFactors = FALSE, check.names = FALSE)
}

safe_adj_r2 <- function(mod) {
  r2 <- tryCatch(vegan::RsquareAdj(mod)$adj.r.squared, error = function(e) NA_real_)
  if (!is.finite(r2)) return(0)
  max(0, r2)
}

# Two-predictor phase varpart (uses existing helper)
# Three-predictor: pure A|B+C, pure B|A+C, pure C|A+B, residual shared, totals
calc_phase_varpart3 <- function(joined_df,
                                guild_cols,
                                preds = c("NAO_Median_Value", "estimate", "tephra"),
                                phase_breaks,
                                min_n = 5) {
  stopifnot(length(preds) == 3)
  A <- preds[1]; B <- preds[2]; C <- preds[3]

  one_phase <- function(bounds, phase_name) {
    lo <- bounds[1]; hi <- bounds[2]
    sub <- joined_df %>%
      dplyr::filter(age_ce >= lo, age_ce < hi) %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(c(guild_cols, A, B, C)), ~ !is.na(.x)))

    comps <- c(
      paste0("Pure ", A), paste0("Pure ", B), paste0("Pure ", C),
      "Shared", "Total_AB", "Total_ABC"
    )
    if (nrow(sub) < min_n) {
      return(tibble::tibble(
        phase = phase_name, component = comps,
        value = rep(NA_real_, length(comps)), n_bins = nrow(sub)
      ))
    }

    # Tephra must vary within the phase slice
    if (length(unique(sub[[C]])) < 2) {
      return(tibble::tibble(
        phase = phase_name, component = comps,
        value = rep(NA_real_, length(comps)), n_bins = nrow(sub)
      ))
    }

    comm <- as.matrix(sub[, guild_cols])
    envA <- as.data.frame(sub[, A, drop = FALSE])
    envB <- as.data.frame(sub[, B, drop = FALSE])
    envC <- as.data.frame(sub[, C, drop = FALSE])
    envAB <- as.data.frame(sub[, c(A, B), drop = FALSE])
    envBC <- as.data.frame(sub[, c(B, C), drop = FALSE])
    envAC <- as.data.frame(sub[, c(A, C), drop = FALSE])
    envABC <- as.data.frame(sub[, c(A, B, C), drop = FALSE])

    pure_A <- safe_adj_r2(vegan::rda(comm, envA, envBC))
    pure_B <- safe_adj_r2(vegan::rda(comm, envB, envAC))
    pure_C <- safe_adj_r2(vegan::rda(comm, envC, envAB))
    total_ab <- safe_adj_r2(vegan::rda(comm, envAB))
    total_abc <- safe_adj_r2(vegan::rda(comm, envABC))
    shared <- max(0, total_abc - pure_A - pure_B - pure_C)

    tibble::tibble(
      phase = phase_name,
      component = comps,
      value = c(pure_A, pure_B, pure_C, shared, total_ab, total_abc),
      n_bins = nrow(sub)
    )
  }

  purrr::imap_dfr(phase_breaks, one_phase)
}

tephra_all <- read_tephra() %>%
  dplyr::mutate(
    lake = revision_normalize_lake_names(lake),
    age_ce = as.numeric(age_ce)
  )

# Defensive: drop any residual regional rows
tephra_all <- tephra_all %>%
  dplyr::filter(!grepl("^none", lake, ignore.case = TRUE))

if ("include_sensitivity" %in% names(tephra_all)) {
  tephra_all$include_sensitivity <- as.logical(tephra_all$include_sensitivity)
} else if ("status" %in% names(tephra_all)) {
  tephra_all$include_sensitivity <-
    grepl("^(SUPPORTED|TENTATIVE)", tephra_all$status, ignore.case = TRUE) &
    is.finite(tephra_all$age_ce)
} else {
  tephra_all$include_sensitivity <- is.finite(tephra_all$age_ce)
}

# Lake-specific SUPPORTED/TENTATIVE tephras only
tephra <- tephra_all %>%
  dplyr::filter(include_sensitivity, is.finite(age_ce))

message("Tephra rows (lake-specific inventory): ", nrow(tephra_all),
        "; used as predictor ages: ", nrow(tephra),
        " across ", length(unique(tephra$lake)), " lakes")

bin_size <- 30
tephra_bins <- tephra %>%
  dplyr::mutate(
    age_bin = floor(age_ce / bin_size) * bin_size,
    influence_lo = age_bin - bin_size,
    influence_hi = age_bin + bin_size
  )
revision_write_csv(tephra_bins, "outputs/revision/A6_tephra_bins.csv")

# Regional 30-yr bins influenced by any retained lake-specific tephra (±1 bin)
influence_ages <- sort(unique(c(
  tephra_bins$age_bin,
  tephra_bins$influence_lo,
  tephra_bins$influence_hi
)))
influence_ages <- influence_ages[is.finite(influence_ages)]

# Intensity: count of unique lakes with a tephra whose ±1-bin window covers the bin
tephra_intensity_by_age <- tephra_bins %>%
  dplyr::rowwise() %>%
  dplyr::summarise(
    lake = lake,
    age_ce = list(c(age_bin, influence_lo, influence_hi)),
    .groups = "drop"
  ) %>%
  tidyr::unnest(age_ce) %>%
  dplyr::distinct(age_ce, lake) %>%
  dplyr::count(age_ce, name = "tephra_n_lakes")

joined <- shared$joined_df_30yr %>%
  dplyr::left_join(tephra_intensity_by_age, by = "age_ce") %>%
  dplyr::mutate(
    tephra = as.integer(age_ce %in% influence_ages),
    tephra_n_lakes = dplyr::coalesce(tephra_n_lakes, 0L)
  )

flagged <- joined %>% dplyr::filter(tephra == 1)
revision_write_csv(
  flagged %>% dplyr::select(age_ce, n_lakes, tephra, tephra_n_lakes),
  "outputs/revision/A6_flagged_bins.csv"
)

guild_cols <- revision_guild_cols
env_ab <- c("NAO_Median_Value", "estimate")
phase_breaks <- revision_phase_breaks

# ---- Primary: with vs without tephra predictor ----
vp_base <- calc_phase_varpart(
  joined_df = joined,
  guild_cols = guild_cols,
  env_predictors = env_ab,
  phase_breaks = phase_breaks,
  min_n = 5
) %>% dplyr::mutate(scenario = "NAO_veg")

vp_tep <- calc_phase_varpart3(
  joined_df = joined,
  guild_cols = guild_cols,
  preds = c("NAO_Median_Value", "estimate", "tephra"),
  phase_breaks = phase_breaks,
  min_n = 5
) %>% dplyr::mutate(scenario = "NAO_veg_tephra")

# Optional intensity alternative (count of lakes with tephra influence)
vp_int <- calc_phase_varpart3(
  joined_df = joined,
  guild_cols = guild_cols,
  preds = c("NAO_Median_Value", "estimate", "tephra_n_lakes"),
  phase_breaks = phase_breaks,
  min_n = 5
) %>%
  dplyr::mutate(
    scenario = "NAO_veg_tephra_intensity",
    component = dplyr::recode(
      component,
      "Pure tephra_n_lakes" = "Pure tephra"
    )
  )

# Optional cheap exclusion diagnostic (not the main result)
joined_ex <- joined %>% dplyr::filter(tephra == 0)
vp_ex <- calc_phase_varpart(
  joined_df = joined_ex,
  guild_cols = guild_cols,
  env_predictors = env_ab,
  phase_breaks = phase_breaks,
  min_n = 5
) %>% dplyr::mutate(scenario = "exclusion_diagnostic")

vp <- dplyr::bind_rows(vp_base, vp_tep, vp_int, vp_ex)
revision_write_csv(vp, "outputs/revision/A6_varpart_tephra_sensitivity.csv")

# Delta: change in pure NAO / vegetation when tephra is added as predictor
base_wide <- vp_base %>%
  dplyr::select(phase, component, value_base = value, n_bins_base = n_bins)
tep_wide <- vp_tep %>%
  dplyr::select(phase, component, value_tep = value, n_bins_tep = n_bins)

delta <- dplyr::full_join(base_wide, tep_wide, by = c("phase", "component")) %>%
  dplyr::mutate(
    delta = value_tep - value_base,
    n_bins = dplyr::coalesce(n_bins_tep, n_bins_base)
  ) %>%
  dplyr::select(phase, component, value_base, value_tep, delta, n_bins)
revision_write_csv(delta, "outputs/revision/A6_varpart_delta.csv")
print(as.data.frame(delta), row.names = FALSE)

# Annotation table for Fig. 2 overlay
src_col <- if ("primary_ref" %in% names(tephra_bins)) {
  ifelse(!is.na(tephra_bins$primary_ref) & tephra_bins$primary_ref != "",
         as.character(tephra_bins$primary_ref),
         as.character(tephra_bins$source))
} else {
  as.character(tephra_bins$source)
}
annot <- tephra_bins %>%
  dplyr::mutate(source = src_col) %>%
  dplyr::distinct(age_ce, event, source) %>%
  dplyr::arrange(age_ce)
revision_write_csv(annot, "outputs/revision/A6_fig2_eruption_annotations.csv")

# Main figure: baseline vs tephra-as-predictor (pure fractions only)
plot_comps <- c("Pure NAO_Median_Value", "Pure estimate", "Pure tephra", "Shared")
vp_plot <- vp %>%
  dplyr::filter(scenario %in% c("NAO_veg", "NAO_veg_tephra")) %>%
  dplyr::filter(component %in% plot_comps) %>%
  dplyr::filter(!is.na(value)) %>%
  dplyr::mutate(
    scenario = factor(
      scenario,
      levels = c("NAO_veg", "NAO_veg_tephra"),
      labels = c("NAO + vegetation", "NAO + vegetation + tephra")
    ),
    component = dplyr::recode(
      component,
      "Pure NAO_Median_Value" = "Pure NAO",
      "Pure estimate" = "Pure vegetation"
    )
  )

p <- ggplot(vp_plot, aes(x = phase, y = value, fill = component)) +
  geom_col(position = "stack") +
  facet_wrap(~scenario) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "A6: Varpart with tephra as predictor (±1×30-yr bin indicator)",
    y = "Adj. R²",
    fill = NULL
  )
ggsave("outputs/revision/figures/A6_tephra_varpart.png", p, width = 9, height = 5, dpi = 150)

message(
  "A6 complete (", nrow(tephra), " tephra ages; ",
  length(influence_ages), " influenced 30-yr bin centres; ",
  nrow(flagged), "/", nrow(joined), " community bins with tephra=1; ",
  "mean tephra prevalence = ", round(mean(joined$tephra), 3), ")"
)
