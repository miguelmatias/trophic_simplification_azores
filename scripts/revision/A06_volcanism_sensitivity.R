# A6 — Volcanism sensitivity: lake-specific tephra as a predictor
# Alternative Figure 5 layout (2 columns × 3 rows):
#   Left:  lake-level NAO + vegetation with Condition(lake)
#   Right: same + lake-specific tephra predictor with Condition(lake)
# Panels: (a) historical phases, (b) moving window, (c) Climate vs VegChange
# Tephra is coded per lake×bin (1 if within ±1×30-yr of that lake's
# SUPPORTED/TENTATIVE tephra ages); no regional tephra curve.
# Significance: permutation tests + faded NS fractions (Fig. 5 style).
source("scripts/revision/_bootstrap.R")
shared <- revision_bootstrap()

suppressPackageStartupMessages({
  library(patchwork)
  library(ggplot2)
})

message("=== A6 alt Fig. 5: lake-level standard vs +tephra (Condition(lake)) ===")

read_tephra <- function(path = "data/revision/tephra_events.csv") {
  if (requireNamespace("readr", quietly = TRUE)) {
    return(as.data.frame(
      readr::read_csv(path, show_col_types = FALSE,
                      locale = readr::locale(encoding = "UTF-8")),
      stringsAsFactors = FALSE
    ))
  }
  raw <- readBin(path, what = "raw", n = file.info(path)$size)
  if (length(raw) >= 3 && raw[1] == as.raw(0xef) &&
      raw[2] == as.raw(0xbb) && raw[3] == as.raw(0xbf)) {
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

# ---- Tephra inventory (lake-specific only) ----
tephra_all <- read_tephra() %>%
  dplyr::mutate(
    lake = revision_normalize_lake_names(lake),
    age_ce = as.numeric(age_ce)
  ) %>%
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

tephra <- tephra_all %>%
  dplyr::filter(include_sensitivity, is.finite(age_ce))

message(
  "Tephra rows (lake-specific inventory): ", nrow(tephra_all),
  "; used as predictor ages: ", nrow(tephra),
  " across ", length(unique(tephra$lake)), " lakes"
)

bin_size <- 30
tephra_bins <- tephra %>%
  dplyr::mutate(
    age_bin = floor(age_ce / bin_size) * bin_size,
    influence_lo = age_bin - bin_size,
    influence_hi = age_bin + bin_size
  )
revision_write_csv(tephra_bins, "outputs/revision/A6_tephra_bins.csv")

# Expanded lake × influenced-bin lookup (centre ±1 bin)
tephra_lake_ages <- tephra_bins %>%
  dplyr::rowwise() %>%
  dplyr::summarise(
    lake = lake,
    age_ce = list(c(age_bin, influence_lo, influence_hi)),
    .groups = "drop"
  ) %>%
  tidyr::unnest(age_ce) %>%
  dplyr::distinct(lake, age_ce) %>%
  dplyr::mutate(tephra = 1L)

revision_write_csv(tephra_lake_ages, "outputs/revision/A6_tephra_lake_bin_matrix.csv")

# ---- Lake-level community + regional env + lake-specific tephra ----
guild_cols <- revision_guild_cols
phase_breaks <- revision_phase_breaks
A <- "NAO_Median_Value"
B <- "estimate"
C <- "tephra"

lake_env <- shared$bio_lake_30 %>%
  dplyr::left_join(shared$nao_30, by = "age_ce") %>%
  dplyr::left_join(shared$veg_30 %>% dplyr::select(age_ce, estimate), by = "age_ce") %>%
  dplyr::left_join(tephra_lake_ages, by = c("lake", "age_ce")) %>%
  dplyr::mutate(
    tephra = dplyr::coalesce(tephra, 0L),
    lake = factor(lake)
  )

flagged <- lake_env %>% dplyr::filter(tephra == 1L)
revision_write_csv(
  flagged %>% dplyr::select(lake, age_ce, tephra, dplyr::all_of(c(A, B))),
  "outputs/revision/A6_flagged_bins.csv"
)

message(
  "Lake×bin rows: ", nrow(lake_env),
  "; tephra=1: ", nrow(flagged),
  " (", round(100 * mean(lake_env$tephra), 1), "%); ",
  "lakes with any tephra flag: ",
  length(unique(flagged$lake))
)

# Minimum rows per phase for lake-level models (A3 used 8)
min_n <- 8L
# Published Fig. 5 used 9999; revision A-scripts use 999 for runtime
n_perm <- as.integer(Sys.getenv("A6_N_PERM", unset = "999"))
set.seed(1)
message("n_perm = ", n_perm, "; min_n = ", min_n)

# ---- Phase varpart helpers with Condition(lake) ----
calc_phase_varpart_lake2 <- function(df, phase_breaks, min_n = 8) {
  purrr::imap_dfr(phase_breaks, function(bounds, phase_name) {
    lo <- bounds[1]; hi <- bounds[2]
    sub <- df %>%
      dplyr::filter(age_ce >= lo, age_ce < hi) %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(c(guild_cols, A, B)), ~ !is.na(.x)))

    comps <- c(paste0("Pure ", A), paste0("Pure ", B), "Shared", "Total_AB")
    if (nrow(sub) < min_n || dplyr::n_distinct(sub$lake) < 2) {
      return(tibble::tibble(
        phase = phase_name, component = comps,
        value = rep(NA_real_, length(comps)),
        n_bins = nrow(sub), n_lakes = dplyr::n_distinct(sub$lake),
        tephra_prevalence = mean(sub$tephra, na.rm = TRUE)
      ))
    }

    comm <- as.matrix(sub[, guild_cols])
    envA <- data.frame(NAO_Median_Value = sub[[A]])
    envB <- data.frame(estimate = sub[[B]])
    lake_df <- data.frame(lake = factor(sub$lake))

    pure_A <- safe_adj_r2(vegan::rda(comm, envA, cbind(envB, lake_df)))
    pure_B <- safe_adj_r2(vegan::rda(comm, envB, cbind(envA, lake_df)))
    total_ab <- safe_adj_r2(vegan::rda(comm, cbind(envA, envB), lake_df))
    shared <- max(0, total_ab - pure_A - pure_B)

    tibble::tibble(
      phase = phase_name,
      component = comps,
      value = c(pure_A, pure_B, shared, total_ab),
      n_bins = nrow(sub),
      n_lakes = dplyr::n_distinct(sub$lake),
      tephra_prevalence = mean(sub$tephra, na.rm = TRUE)
    )
  })
}

calc_phase_varpart_lake3 <- function(df, phase_breaks, min_n = 8) {
  purrr::imap_dfr(phase_breaks, function(bounds, phase_name) {
    lo <- bounds[1]; hi <- bounds[2]
    sub <- df %>%
      dplyr::filter(age_ce >= lo, age_ce < hi) %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(c(guild_cols, A, B, C)), ~ !is.na(.x)))

    comps <- c(
      paste0("Pure ", A), paste0("Pure ", B), paste0("Pure ", C),
      "Shared", "Total_AB", "Total_ABC"
    )
    n_lakes <- dplyr::n_distinct(sub$lake)
    tep_prev <- mean(sub$tephra, na.rm = TRUE)

    if (nrow(sub) < min_n || n_lakes < 2 || length(unique(sub[[C]])) < 2) {
      return(tibble::tibble(
        phase = phase_name, component = comps,
        value = rep(NA_real_, length(comps)),
        n_bins = nrow(sub), n_lakes = n_lakes,
        tephra_prevalence = tep_prev
      ))
    }

    comm <- as.matrix(sub[, guild_cols])
    envA <- data.frame(NAO_Median_Value = sub[[A]])
    envB <- data.frame(estimate = sub[[B]])
    envC <- data.frame(tephra = sub[[C]])
    lake_df <- data.frame(lake = factor(sub$lake))

    pure_A <- safe_adj_r2(vegan::rda(comm, envA, cbind(envB, envC, lake_df)))
    pure_B <- safe_adj_r2(vegan::rda(comm, envB, cbind(envA, envC, lake_df)))
    pure_C <- safe_adj_r2(vegan::rda(comm, envC, cbind(envA, envB, lake_df)))
    total_ab <- safe_adj_r2(vegan::rda(comm, cbind(envA, envB), lake_df))
    total_abc <- safe_adj_r2(vegan::rda(comm, cbind(envA, envB, envC), lake_df))
    shared <- max(0, total_abc - pure_A - pure_B - pure_C)

    tibble::tibble(
      phase = phase_name,
      component = comps,
      value = c(pure_A, pure_B, pure_C, shared, total_ab, total_abc),
      n_bins = nrow(sub),
      n_lakes = n_lakes,
      tephra_prevalence = tep_prev
    )
  })
}

# ---- Significance (mirror Fig. 5 / A11; Condition(lake) throughout) ----
phase_sig2 <- purrr::imap_dfr(phase_breaks, function(bounds, phase_name) {
  lo <- bounds[1]; hi <- bounds[2]
  sub <- lake_env %>%
    dplyr::filter(age_ce >= lo, age_ce < hi) %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(c(guild_cols, A, B)), ~ !is.na(.x)))

  n_bins <- nrow(sub)
  n_lakes <- dplyr::n_distinct(sub$lake)
  if (n_bins < min_n || n_lakes < 2) {
    return(tibble::tibble(
      phase = phase_name, n_bins = n_bins, n_lakes = n_lakes,
      p_full = NA_real_, p_pure_A = NA_real_, p_pure_B = NA_real_
    ))
  }

  comm <- as.matrix(sub[, guild_cols])
  env <- data.frame(
    NAO_Median_Value = sub[[A]],
    estimate = sub[[B]],
    lake = factor(sub$lake)
  )

  mod_full <- vegan::rda(
    as.formula(paste("comm ~", A, "+", B, "+ Condition(lake)")),
    data = env
  )
  mod_pure_A <- vegan::rda(
    as.formula(paste("comm ~", A, "+ Condition(", B, " + lake)")),
    data = env
  )
  mod_pure_B <- vegan::rda(
    as.formula(paste("comm ~", B, "+ Condition(", A, " + lake)")),
    data = env
  )

  tibble::tibble(
    phase = phase_name, n_bins = n_bins, n_lakes = n_lakes,
    p_full = safe_anova_p(mod_full, permutations = n_perm),
    p_pure_A = safe_anova_p(mod_pure_A, permutations = n_perm),
    p_pure_B = safe_anova_p(mod_pure_B, permutations = n_perm)
  )
})

message("Computing 3-predictor phase significance...")
phase_sig3 <- purrr::imap_dfr(phase_breaks, function(bounds, phase_name) {
  lo <- bounds[1]; hi <- bounds[2]
  sub <- lake_env %>%
    dplyr::filter(age_ce >= lo, age_ce < hi) %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(c(guild_cols, A, B, C)), ~ !is.na(.x)))

  n_bins <- nrow(sub)
  n_lakes <- dplyr::n_distinct(sub$lake)
  tep_varies <- length(unique(sub[[C]])) >= 2

  if (n_bins < min_n || n_lakes < 2 || !tep_varies) {
    return(tibble::tibble(
      phase = phase_name, n_bins = n_bins, n_lakes = n_lakes,
      tephra_varies = tep_varies,
      p_full = NA_real_, p_pure_A = NA_real_, p_pure_B = NA_real_,
      p_pure_C = NA_real_
    ))
  }

  comm <- as.matrix(sub[, guild_cols])
  env <- data.frame(
    NAO_Median_Value = sub[[A]],
    estimate = sub[[B]],
    tephra = sub[[C]],
    lake = factor(sub$lake)
  )

  mod_full <- vegan::rda(
    as.formula(paste("comm ~", A, "+", B, "+", C, "+ Condition(lake)")),
    data = env
  )
  mod_pure_A <- vegan::rda(
    as.formula(paste("comm ~", A, "+ Condition(", B, " + ", C, " + lake)")),
    data = env
  )
  mod_pure_B <- vegan::rda(
    as.formula(paste("comm ~", B, "+ Condition(", A, " + ", C, " + lake)")),
    data = env
  )
  mod_pure_C <- vegan::rda(
    as.formula(paste("comm ~", C, "+ Condition(", A, " + ", B, " + lake)")),
    data = env
  )

  tibble::tibble(
    phase = phase_name, n_bins = n_bins, n_lakes = n_lakes,
    tephra_varies = TRUE,
    p_full = safe_anova_p(mod_full, permutations = n_perm),
    p_pure_A = safe_anova_p(mod_pure_A, permutations = n_perm),
    p_pure_B = safe_anova_p(mod_pure_B, permutations = n_perm),
    p_pure_C = safe_anova_p(mod_pure_C, permutations = n_perm)
  )
})

# ---- Effect sizes ----
message("Computing phase varpart effect sizes...")
vp_base <- calc_phase_varpart_lake2(lake_env, phase_breaks, min_n = min_n) %>%
  dplyr::mutate(scenario = "NAO_veg")
vp_tep <- calc_phase_varpart_lake3(lake_env, phase_breaks, min_n = min_n) %>%
  dplyr::mutate(scenario = "NAO_veg_tephra")

vp <- dplyr::bind_rows(vp_base, vp_tep)
revision_write_csv(vp, "outputs/revision/A6_varpart_tephra_sensitivity.csv")
revision_write_csv(phase_sig2, "outputs/revision/A6_phase_sig_base.csv")
revision_write_csv(phase_sig3, "outputs/revision/A6_phase_sig_tephra.csv")

# Delta table
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

# SI-ready results with p-values
sig_join <- phase_sig3 %>%
  dplyr::select(phase, p_full, p_pure_A, p_pure_B, p_pure_C, tephra_varies)
effects_wide <- vp_tep %>%
  dplyr::filter(component %in% c(
    paste0("Pure ", A), paste0("Pure ", B), paste0("Pure ", C),
    "Shared", "Total_ABC"
  )) %>%
  dplyr::select(phase, component, value, n_bins, n_lakes, tephra_prevalence) %>%
  tidyr::pivot_wider(names_from = component, values_from = value)

results_table <- effects_wide %>%
  dplyr::left_join(sig_join, by = "phase") %>%
  dplyr::left_join(
    phase_sig2 %>%
      dplyr::select(phase, p_full_base = p_full,
                    p_pure_A_base = p_pure_A, p_pure_B_base = p_pure_B),
    by = "phase"
  )
revision_write_csv(results_table, "outputs/revision/A6_varpart_results_with_p.csv")
print(as.data.frame(results_table), row.names = FALSE)

# Annotation table for Fig. 2 overlay (ages only; not a predictor curve)
src_col <- if ("primary_ref" %in% names(tephra_bins)) {
  ifelse(!is.na(tephra_bins$primary_ref) & tephra_bins$primary_ref != "",
         as.character(tephra_bins$primary_ref),
         as.character(tephra_bins$source))
} else {
  as.character(tephra_bins$source)
}
annot <- tephra_bins %>%
  dplyr::mutate(source = src_col) %>%
  dplyr::distinct(lake, age_ce, event, source) %>%
  dplyr::arrange(lake, age_ce)
revision_write_csv(annot, "outputs/revision/A6_fig2_eruption_annotations.csv")

# =====================================================================
# Moving-window lake-level varpart (both columns Condition(lake))
# =====================================================================
window_size <- 300
step_size   <- 30
min_year <- floor(min(lake_env$age_ce, na.rm = TRUE))
max_year <- ceiling(max(lake_env$age_ce, na.rm = TRUE))
window_starts <- seq(min_year, max_year - window_size, by = step_size)
message(
  "Moving windows: ", length(window_starts),
  " (size = ", window_size, ", step = ", step_size, ")"
)

calc_window_varpart_lake2 <- function(df, window_start, window_size, min_n = 8) {
  window_end <- window_start + window_size
  sub <- df %>%
    dplyr::filter(age_ce >= window_start, age_ce < window_end) %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(c(guild_cols, A, B)), ~ !is.na(.x)))

  mid <- (window_start + window_end) / 2
  comps <- c("Pure Climate", "Pure Vegetation", "Shared", "Total_AB")
  n_lakes <- dplyr::n_distinct(sub$lake)

  if (nrow(sub) < min_n || n_lakes < 2) {
    return(tibble::tibble(
      window_start = window_start, window_end = window_end, midpoint = mid,
      component = comps, value = rep(NA_real_, length(comps)),
      n_bins = nrow(sub), n_lakes = n_lakes,
      tephra_prevalence = mean(sub$tephra, na.rm = TRUE)
    ))
  }

  comm <- as.matrix(sub[, guild_cols])
  envA <- data.frame(NAO_Median_Value = sub[[A]])
  envB <- data.frame(estimate = sub[[B]])
  lake_df <- data.frame(lake = factor(sub$lake))

  pure_A <- safe_adj_r2(vegan::rda(comm, envA, cbind(envB, lake_df)))
  pure_B <- safe_adj_r2(vegan::rda(comm, envB, cbind(envA, lake_df)))
  total_ab <- safe_adj_r2(vegan::rda(comm, cbind(envA, envB), lake_df))
  shared <- max(0, total_ab - pure_A - pure_B)

  tibble::tibble(
    window_start = window_start, window_end = window_end, midpoint = mid,
    component = comps,
    value = c(pure_A, pure_B, shared, total_ab),
    n_bins = nrow(sub), n_lakes = n_lakes,
    tephra_prevalence = mean(sub$tephra, na.rm = TRUE)
  )
}

calc_window_varpart_lake3 <- function(df, window_start, window_size, min_n = 8) {
  window_end <- window_start + window_size
  sub <- df %>%
    dplyr::filter(age_ce >= window_start, age_ce < window_end) %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(c(guild_cols, A, B, C)), ~ !is.na(.x)))

  mid <- (window_start + window_end) / 2
  comps <- c(
    "Pure Climate", "Pure Vegetation", "Pure tephra",
    "Shared", "Total_AB", "Total_ABC"
  )
  n_lakes <- dplyr::n_distinct(sub$lake)
  tep_prev <- mean(sub$tephra, na.rm = TRUE)
  tep_varies <- length(unique(sub[[C]])) >= 2

  if (nrow(sub) < min_n || n_lakes < 2 || !tep_varies) {
    return(tibble::tibble(
      window_start = window_start, window_end = window_end, midpoint = mid,
      component = comps, value = rep(NA_real_, length(comps)),
      n_bins = nrow(sub), n_lakes = n_lakes,
      tephra_prevalence = tep_prev
    ))
  }

  comm <- as.matrix(sub[, guild_cols])
  envA <- data.frame(NAO_Median_Value = sub[[A]])
  envB <- data.frame(estimate = sub[[B]])
  envC <- data.frame(tephra = sub[[C]])
  lake_df <- data.frame(lake = factor(sub$lake))

  pure_A <- safe_adj_r2(vegan::rda(comm, envA, cbind(envB, envC, lake_df)))
  pure_B <- safe_adj_r2(vegan::rda(comm, envB, cbind(envA, envC, lake_df)))
  pure_C <- safe_adj_r2(vegan::rda(comm, envC, cbind(envA, envB, lake_df)))
  total_ab <- safe_adj_r2(vegan::rda(comm, cbind(envA, envB), lake_df))
  total_abc <- safe_adj_r2(vegan::rda(comm, cbind(envA, envB, envC), lake_df))
  shared <- max(0, total_abc - pure_A - pure_B - pure_C)

  tibble::tibble(
    window_start = window_start, window_end = window_end, midpoint = mid,
    component = comps,
    value = c(pure_A, pure_B, pure_C, shared, total_ab, total_abc),
    n_bins = nrow(sub), n_lakes = n_lakes,
    tephra_prevalence = tep_prev
  )
}

calc_window_sig_lake2 <- function(df, window_start, window_size, min_n = 8, n_perm = 999) {
  window_end <- window_start + window_size
  sub <- df %>%
    dplyr::filter(age_ce >= window_start, age_ce < window_end) %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(c(guild_cols, A, B)), ~ !is.na(.x)))

  mid <- (window_start + window_end) / 2
  n_bins <- nrow(sub)
  n_lakes <- dplyr::n_distinct(sub$lake)

  if (n_bins < min_n || n_lakes < 2) {
    return(tibble::tibble(
      window_start = window_start, window_end = window_end, midpoint = mid,
      n_bins = n_bins, n_lakes = n_lakes,
      p_full = NA_real_, p_pure_A = NA_real_, p_pure_B = NA_real_
    ))
  }

  comm <- as.matrix(sub[, guild_cols])
  env <- data.frame(
    NAO_Median_Value = sub[[A]],
    estimate = sub[[B]],
    lake = factor(sub$lake)
  )

  mod_full <- vegan::rda(
    as.formula(paste("comm ~", A, "+", B, "+ Condition(lake)")), data = env
  )
  mod_pure_A <- vegan::rda(
    as.formula(paste("comm ~", A, "+ Condition(", B, " + lake)")), data = env
  )
  mod_pure_B <- vegan::rda(
    as.formula(paste("comm ~", B, "+ Condition(", A, " + lake)")), data = env
  )

  tibble::tibble(
    window_start = window_start, window_end = window_end, midpoint = mid,
    n_bins = n_bins, n_lakes = n_lakes,
    p_full = safe_anova_p(mod_full, permutations = n_perm),
    p_pure_A = safe_anova_p(mod_pure_A, permutations = n_perm),
    p_pure_B = safe_anova_p(mod_pure_B, permutations = n_perm)
  )
}

calc_window_sig_lake3 <- function(df, window_start, window_size, min_n = 8, n_perm = 999) {
  window_end <- window_start + window_size
  sub <- df %>%
    dplyr::filter(age_ce >= window_start, age_ce < window_end) %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(c(guild_cols, A, B, C)), ~ !is.na(.x)))

  mid <- (window_start + window_end) / 2
  n_bins <- nrow(sub)
  n_lakes <- dplyr::n_distinct(sub$lake)
  tep_varies <- length(unique(sub[[C]])) >= 2

  if (n_bins < min_n || n_lakes < 2 || !tep_varies) {
    return(tibble::tibble(
      window_start = window_start, window_end = window_end, midpoint = mid,
      n_bins = n_bins, n_lakes = n_lakes, tephra_varies = tep_varies,
      p_full = NA_real_, p_pure_A = NA_real_, p_pure_B = NA_real_,
      p_pure_C = NA_real_
    ))
  }

  comm <- as.matrix(sub[, guild_cols])
  env <- data.frame(
    NAO_Median_Value = sub[[A]],
    estimate = sub[[B]],
    tephra = sub[[C]],
    lake = factor(sub$lake)
  )

  mod_full <- vegan::rda(
    as.formula(paste("comm ~", A, "+", B, "+", C, "+ Condition(lake)")), data = env
  )
  mod_pure_A <- vegan::rda(
    as.formula(paste("comm ~", A, "+ Condition(", B, " + ", C, " + lake)")), data = env
  )
  mod_pure_B <- vegan::rda(
    as.formula(paste("comm ~", B, "+ Condition(", A, " + ", C, " + lake)")), data = env
  )
  mod_pure_C <- vegan::rda(
    as.formula(paste("comm ~", C, "+ Condition(", A, " + ", B, " + lake)")), data = env
  )

  tibble::tibble(
    window_start = window_start, window_end = window_end, midpoint = mid,
    n_bins = n_bins, n_lakes = n_lakes, tephra_varies = TRUE,
    p_full = safe_anova_p(mod_full, permutations = n_perm),
    p_pure_A = safe_anova_p(mod_pure_A, permutations = n_perm),
    p_pure_B = safe_anova_p(mod_pure_B, permutations = n_perm),
    p_pure_C = safe_anova_p(mod_pure_C, permutations = n_perm)
  )
}

message("Computing moving-window effect sizes (base + tephra)...")
win_base <- purrr::map_dfr(
  window_starts,
  ~ calc_window_varpart_lake2(lake_env, .x, window_size, min_n = min_n)
) %>% dplyr::mutate(scenario = "NAO_veg")

win_tep <- purrr::map_dfr(
  window_starts,
  ~ calc_window_varpart_lake3(lake_env, .x, window_size, min_n = min_n)
) %>% dplyr::mutate(scenario = "NAO_veg_tephra")

win_vp <- dplyr::bind_rows(win_base, win_tep)
revision_write_csv(win_vp, "outputs/revision/A6_window_varpart.csv")

message("Computing moving-window significance (base; n_perm = ", n_perm, ")...")
win_sig2 <- purrr::map_dfr(
  window_starts,
  ~ calc_window_sig_lake2(lake_env, .x, window_size, min_n = min_n, n_perm = n_perm)
) %>% dplyr::mutate(scenario = "NAO_veg")

message("Computing moving-window significance (+tephra; n_perm = ", n_perm, ")...")
win_sig3 <- purrr::map_dfr(
  window_starts,
  ~ calc_window_sig_lake3(lake_env, .x, window_size, min_n = min_n, n_perm = n_perm)
) %>% dplyr::mutate(scenario = "NAO_veg_tephra")

win_sig <- dplyr::bind_rows(
  win_sig2 %>% dplyr::mutate(p_pure_C = NA_real_, tephra_varies = NA),
  win_sig3
)
revision_write_csv(win_sig, "outputs/revision/A6_window_sig.csv")

# ---- Plot data: phases ----
alpha_sig <- 1.00
alpha_ns  <- 0.25

component_colors <- c(
  "Pure Climate" = "#F4B942",
  "Pure Vegetation" = "#97D8C4",
  "Pure tephra" = "#8FBC8F",
  "Shared" = "#58B0E6"
)
effect_colors_split <- c(
  "VegChange > Climate (pre-750)"  = "#054A29",
  "VegChange > Climate (post-750)" = "#6DC6B6",
  "Climate > VegChange"            = "#F4B942"
)

phase_flags <- dplyr::bind_rows(
  phase_sig2 %>%
    dplyr::transmute(
      phase, scenario = "NAO_veg",
      sig_pure_A = !is.na(p_pure_A) & p_pure_A < 0.05,
      sig_pure_B = !is.na(p_pure_B) & p_pure_B < 0.05,
      sig_pure_C = FALSE,
      sig_full   = !is.na(p_full) & p_full < 0.05
    ),
  phase_sig3 %>%
    dplyr::transmute(
      phase, scenario = "NAO_veg_tephra",
      sig_pure_A = !is.na(p_pure_A) & p_pure_A < 0.05,
      sig_pure_B = !is.na(p_pure_B) & p_pure_B < 0.05,
      sig_pure_C = !is.na(p_pure_C) & p_pure_C < 0.05,
      sig_full   = !is.na(p_full) & p_full < 0.05
    )
)

phase_plot <- vp %>%
  dplyr::filter(component %in% c(
    paste0("Pure ", A), paste0("Pure ", B), paste0("Pure ", C), "Shared"
  )) %>%
  dplyr::filter(!is.na(value), value > 0) %>%
  dplyr::left_join(phase_flags, by = c("phase", "scenario")) %>%
  dplyr::mutate(
    component_label = dplyr::recode(
      component,
      "Pure NAO_Median_Value" = "Pure Climate",
      "Pure estimate" = "Pure Vegetation",
      "Pure tephra" = "Pure tephra",
      "Shared" = "Shared"
    ),
    # Left column never shows tephra
    keep = !(scenario == "NAO_veg" & component_label == "Pure tephra")
  ) %>%
  dplyr::filter(keep) %>%
  dplyr::mutate(
    alpha_flag = dplyr::case_when(
      component_label == "Pure Climate" ~ ifelse(sig_pure_A, alpha_sig, alpha_ns),
      component_label == "Pure Vegetation" ~ ifelse(sig_pure_B, alpha_sig, alpha_ns),
      component_label == "Pure tephra" ~ ifelse(sig_pure_C, alpha_sig, alpha_ns),
      component_label == "Shared" ~ alpha_sig,
      TRUE ~ alpha_sig
    ),
    column = factor(
      scenario,
      levels = c("NAO_veg", "NAO_veg_tephra"),
      labels = c("NAO + vegetation", "NAO + vegetation + tephra")
    ),
    component_label = factor(
      component_label,
      levels = c("Shared", "Pure Vegetation", "Pure tephra", "Pure Climate")
    ),
    phase = factor(phase, levels = paste("Phase", 1:5))
  )

# ---- Plot data: windows ----
win_flags <- win_sig %>%
  dplyr::transmute(
    window_start, window_end, midpoint, scenario,
    sig_pure_A = !is.na(p_pure_A) & p_pure_A < 0.05,
    sig_pure_B = !is.na(p_pure_B) & p_pure_B < 0.05,
    sig_pure_C = !is.na(p_pure_C) & p_pure_C < 0.05
  )

win_plot <- win_vp %>%
  dplyr::filter(component %in% c(
    "Pure Climate", "Pure Vegetation", "Pure tephra", "Shared"
  )) %>%
  dplyr::filter(!is.na(value), value > 0) %>%
  dplyr::left_join(win_flags, by = c("window_start", "window_end", "midpoint", "scenario")) %>%
  dplyr::mutate(
    keep = !(scenario == "NAO_veg" & component == "Pure tephra")
  ) %>%
  dplyr::filter(keep) %>%
  dplyr::mutate(
    alpha_flag = dplyr::case_when(
      component == "Pure Climate" ~ ifelse(sig_pure_A, alpha_sig, alpha_ns),
      component == "Pure Vegetation" ~ ifelse(sig_pure_B, alpha_sig, alpha_ns),
      component == "Pure tephra" ~ ifelse(sig_pure_C, alpha_sig, alpha_ns),
      component == "Shared" ~ alpha_sig,
      TRUE ~ alpha_sig
    ),
    column = factor(
      scenario,
      levels = c("NAO_veg", "NAO_veg_tephra"),
      labels = c("NAO + vegetation", "NAO + vegetation + tephra")
    ),
    component = factor(
      component,
      levels = c("Shared", "Pure Vegetation", "Pure tephra", "Pure Climate")
    )
  )

# ---- Plot data: effect difference (Vegetation − Climate) ----
effect_diff <- win_vp %>%
  dplyr::filter(component %in% c("Pure Vegetation", "Pure Climate")) %>%
  dplyr::select(window_start, window_end, midpoint, scenario, component, value) %>%
  tidyr::pivot_wider(names_from = component, values_from = value) %>%
  dplyr::left_join(
    win_sig %>%
      dplyr::select(window_start, window_end, midpoint, scenario, p_pure_A, p_pure_B),
    by = c("window_start", "window_end", "midpoint", "scenario")
  ) %>%
  dplyr::mutate(
    `Pure Vegetation` = dplyr::coalesce(`Pure Vegetation`, 0),
    `Pure Climate` = dplyr::coalesce(`Pure Climate`, 0),
    effect_diff = `Pure Vegetation` - `Pure Climate`,
    effect_group = dplyr::case_when(
      effect_diff > 0 & midpoint < 750  ~ "VegChange > Climate (pre-750)",
      effect_diff > 0 & midpoint >= 750 ~ "VegChange > Climate (post-750)",
      effect_diff < 0                   ~ "Climate > VegChange",
      TRUE                              ~ NA_character_
    ),
    sig_clim = !is.na(p_pure_A) & p_pure_A < 0.05,
    sig_veg  = !is.na(p_pure_B) & p_pure_B < 0.05,
    direction = dplyr::case_when(
      effect_diff > 0 ~ "veg",
      effect_diff < 0 ~ "clim",
      TRUE            ~ NA_character_
    ),
    alpha_flag = dplyr::case_when(
      direction == "veg"  & sig_veg  ~ alpha_sig,
      direction == "veg"  & !sig_veg ~ alpha_ns,
      direction == "clim" & sig_clim ~ alpha_sig,
      direction == "clim" & !sig_clim ~ alpha_ns,
      TRUE ~ alpha_ns
    ),
    column = factor(
      scenario,
      levels = c("NAO_veg", "NAO_veg_tephra"),
      labels = c("NAO + vegetation", "NAO + vegetation + tephra")
    )
  ) %>%
  dplyr::filter(!is.na(effect_group))

revision_write_csv(phase_plot, "outputs/revision/A6_alt_fig5_phase_plot.csv")
revision_write_csv(win_plot, "outputs/revision/A6_alt_fig5_window_plot.csv")
revision_write_csv(effect_diff, "outputs/revision/A6_alt_fig5_effect_diff.csv")

# =====================================================================
# Alternative Figure 5: 2 columns × 3 rows (both lake-conditioned)
# =====================================================================
x_min <- floor(min(lake_env$age_ce, na.rm = TRUE) / 100) * 100
x_max <- ceiling(max(lake_env$age_ce, na.rm = TRUE) / 100) * 100
common_breaks <- seq(x_min, x_max, by = 200)
common_minor  <- seq(x_min, x_max, by = 100)

base_size       <- 9
title_size      <- 10
axis_title_size <- 9
axis_text_size  <- 8
legend_title_sz <- 8
legend_text_sz  <- 7

panel_theme <- ggplot2::theme_minimal(base_size = base_size) +
  ggplot2::theme(
    plot.title        = ggplot2::element_text(size = title_size, face = "bold", hjust = 0),
    axis.title        = ggplot2::element_text(size = axis_title_size),
    axis.text         = ggplot2::element_text(size = axis_text_size),
    strip.text        = ggplot2::element_text(size = title_size, face = "bold"),
    legend.position   = "bottom",
    legend.title      = ggplot2::element_text(size = legend_title_sz),
    legend.text       = ggplot2::element_text(size = legend_text_sz),
    legend.key.height = grid::unit(2.5, "mm"),
    legend.key.width  = grid::unit(5.0, "mm"),
    plot.margin       = ggplot2::margin(2, 4, 2, 2)
  )

# Shared y-limits per row for fair left/right comparison
ylim_a <- c(0, max(phase_plot$value, na.rm = TRUE) * 1.15)
ylim_b <- c(0, max(win_plot$value, na.rm = TRUE) * 1.10)
ylim_c <- range(effect_diff$effect_diff, na.rm = TRUE)
ylim_c <- c(min(ylim_c[1] * 1.1, -0.02), max(ylim_c[2] * 1.1, 0.02))

label_phase <- phase_plot %>% dplyr::filter(value > 0.015)

# Right-column tephra presence ticks (not full-height n/a grey bars).
# Phase ticks: phases with any lake×bin tephra = 1 (prevalence > 0).
# Continuous ticks: unique 30-yr bin centres with tephra = 1 in the analysis matrix.
tephra_col_lvl <- factor(
  "NAO + vegetation + tephra",
  levels = c("NAO + vegetation", "NAO + vegetation + tephra")
)

tep_phase_ticks <- results_table %>%
  dplyr::filter(is.finite(tephra_prevalence), tephra_prevalence > 0) %>%
  dplyr::mutate(
    column = tephra_col_lvl,
    phase = factor(phase, levels = paste("Phase", 1:5)),
    y = ylim_a[2]
  )

# Light text only for non-estimable Phase 5 (no tall grey bar)
na_phase_ann <- phase_sig3 %>%
  dplyr::filter(!tephra_varies) %>%
  dplyr::mutate(
    column = tephra_col_lvl,
    phase = factor(phase, levels = paste("Phase", 1:5)),
    y = ylim_a[2] * 0.55,
    label = "n/a — no\ntephra variation"
  )

tep_age_ticks <- flagged %>%
  dplyr::distinct(age_ce) %>%
  dplyr::arrange(age_ce) %>%
  dplyr::mutate(column = tephra_col_lvl)

revision_write_csv(
  tep_age_ticks %>% dplyr::select(age_ce),
  "outputs/revision/A6_tephra_presence_ages.csv"
)

p_a <- ggplot2::ggplot(
  phase_plot,
  ggplot2::aes(x = phase, y = value, fill = component_label, alpha = alpha_flag)
) +
  ggplot2::geom_col(width = 0.65, position = "stack") +
  ggplot2::geom_text(
    data = label_phase,
    ggplot2::aes(label = paste0(round(value * 100, 1), "%")),
    position = ggplot2::position_stack(vjust = 0.5),
    color = "black", size = 2.1, alpha = 1
  ) +
  ggplot2::geom_point(
    data = tep_phase_ticks,
    ggplot2::aes(x = phase, y = y),
    inherit.aes = FALSE,
    shape = 18, size = 2.4, color = "#2F6B4F"
  ) +
  ggplot2::geom_text(
    data = na_phase_ann,
    ggplot2::aes(x = phase, y = y, label = label),
    inherit.aes = FALSE,
    color = "grey35", size = 2.2, lineheight = 0.9, fontface = "italic"
  ) +
  ggplot2::facet_wrap(~column, nrow = 1) +
  ggplot2::scale_alpha_identity() +
  ggplot2::scale_fill_manual(
    name = "Variance Component",
    values = component_colors,
    breaks = c("Pure Climate", "Pure Vegetation", "Pure tephra", "Shared"),
    drop = FALSE
  ) +
  ggplot2::scale_x_discrete(drop = FALSE) +
  ggplot2::coord_cartesian(ylim = ylim_a, clip = "off") +
  ggplot2::labs(
    x = NULL,
    y = expression("Adjusted "*R^2),
    title = "(a) Historical phases"
  ) +
  panel_theme +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
    plot.margin = ggplot2::margin(6, 4, 2, 2)
  )

p_b <- ggplot2::ggplot(
  win_plot,
  ggplot2::aes(x = midpoint, y = value, fill = component, alpha = alpha_flag)
) +
  ggplot2::geom_col(width = step_size, position = "stack") +
  ggplot2::geom_point(
    data = tep_age_ticks,
    ggplot2::aes(x = age_ce, y = ylim_b[2]),
    inherit.aes = FALSE,
    shape = "|", size = 2.8, color = "#2F6B4F"
  ) +
  ggplot2::facet_wrap(~column, nrow = 1) +
  ggplot2::scale_alpha_identity() +
  ggplot2::scale_fill_manual(
    name = "Variance Component",
    values = component_colors,
    breaks = c("Pure Climate", "Pure Vegetation", "Pure tephra", "Shared"),
    drop = FALSE
  ) +
  ggplot2::scale_x_continuous(
    limits = c(x_min, x_max),
    breaks = common_breaks,
    minor_breaks = common_minor,
    expand = c(0, 0)
  ) +
  ggplot2::coord_cartesian(ylim = ylim_b, clip = "off") +
  ggplot2::labs(
    x = "Year (CE)",
    y = expression("Adjusted "*R^2),
    title = "(b) Moving window"
  ) +
  panel_theme +
  ggplot2::theme(plot.margin = ggplot2::margin(6, 4, 2, 2))

p_c <- ggplot2::ggplot(
  effect_diff,
  ggplot2::aes(x = midpoint, y = effect_diff, fill = effect_group, alpha = alpha_flag)
) +
  ggplot2::geom_col(width = step_size) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
  ggplot2::geom_point(
    data = tep_age_ticks,
    ggplot2::aes(x = age_ce, y = ylim_c[2]),
    inherit.aes = FALSE,
    shape = "|", size = 2.8, color = "#2F6B4F"
  ) +
  ggplot2::facet_wrap(~column, nrow = 1) +
  ggplot2::scale_alpha_identity() +
  ggplot2::scale_fill_manual(
    name = "Dominant Driver",
    values = effect_colors_split,
    breaks = c(
      "Climate > VegChange",
      "VegChange > Climate (post-750)",
      "VegChange > Climate (pre-750)"
    ),
    drop = FALSE
  ) +
  ggplot2::scale_x_continuous(
    limits = c(x_min, x_max),
    breaks = common_breaks,
    minor_breaks = common_minor,
    expand = c(0, 0)
  ) +
  ggplot2::coord_cartesian(ylim = ylim_c, clip = "off") +
  ggplot2::labs(
    x = "Year (CE)",
    y = "Effect size (Vegetation − Climate)",
    title = "(c) Climate vs VegChange",
    caption = paste0(
      "Lake-level alternative Fig. 5 (both columns: Condition(lake)). ",
      "Left = NAO + vegetation; right = + lake-specific tephra ",
      "(tephra = 1 only within ±1×30-yr of that lake's SUPPORTED/TENTATIVE ages; ",
      "no regional tephra curve). Faded = permutation p ≥ 0.05. ",
      "Right-panel ticks mark lake×bins (b–c) or phases (a) with tephra = 1 present; ",
      "Phase 5 n/a text = no tephra variation (model skipped), not a null effect. ",
      "n_perm = ", n_perm, " (published Fig. 5 used 9999)."
    )
  ) +
  panel_theme +
  ggplot2::theme(
    plot.caption = ggplot2::element_text(
      size = 6.5, hjust = 0, colour = "grey30", lineheight = 1.1
    ),
    plot.margin = ggplot2::margin(6, 4, 2, 2)
  )

alt_fig5 <- (p_a / p_b / p_c) +
  patchwork::plot_layout(guides = "collect", heights = c(1, 1.15, 1.15)) &
  ggplot2::theme(legend.position = "bottom")

out_png <- "outputs/revision/figures/A6_alt_figure5_standard_vs_tephra.png"
out_pdf <- "outputs/revision/figures/A6_alt_figure5_standard_vs_tephra.pdf"

ggplot2::ggsave(
  out_png, alt_fig5,
  width = 9.0, height = 10.5, units = "in", dpi = 300, bg = "white"
)
ggplot2::ggsave(
  out_pdf, alt_fig5,
  width = 9.0, height = 10.5, units = "in", device = "pdf", bg = "white"
)

# Keep simpler phase-only extract (same data as row a) for SI continuity
p_phase_only <- p_a +
  ggplot2::labs(
    title = "A6: Lake-level varpart with lake-specific tephra (+ Condition(lake))",
    subtitle = paste0(
      "Faded fractions = non-significant (permutation p ≥ 0.05); ",
      "extract of alt Fig. 5 row (a); diamonds = phases with tephra present"
    )
  )
ggplot2::ggsave(
  "outputs/revision/figures/A6_tephra_varpart.png",
  p_phase_only, width = 9.0, height = 4.2, dpi = 150, bg = "white"
)
ggplot2::ggsave(
  "outputs/revision/figures/A6_tephra_varpart.pdf",
  p_phase_only, width = 9.0, height = 4.2, device = "pdf", bg = "white"
)

message("Saved ", out_png, " and ", out_pdf)
message(
  "A6 complete (", nrow(tephra), " tephra ages; ",
  nrow(flagged), "/", nrow(lake_env), " lake×bins with tephra=1; ",
  "mean tephra prevalence = ", round(mean(lake_env$tephra), 3), ")"
)
