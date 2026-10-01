# A6 — Volcanism sensitivity: lake-specific tephra as a predictor
# Primary test: phase varpart on lake × 30-yr bins with Condition(lake)
#   community ~ NAO + vegetation
#   vs community ~ NAO + vegetation + tephra
# Tephra is coded per lake×bin (1 if within ±1×30-yr of that lake's
# SUPPORTED/TENTATIVE tephra ages); no regional tephra curve.
# Significance: permutation tests + faded NS fractions (Fig. 5 style).
source("scripts/revision/_bootstrap.R")
shared <- revision_bootstrap()

message("=== A6 volcanism sensitivity (lake-level tephra + Condition(lake)) ===")

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

# ---- Figure: baseline vs +tephra with faded NS fractions (Fig. 5 style) ----
alpha_sig <- 1.00
alpha_ns  <- 0.25
shared_rule <- "always"  # match published Fig. 5 phase panel

component_colors <- c(
  "Pure NAO" = "#F4B942",
  "Pure vegetation" = "#97D8C4",
  "Pure tephra" = "#8FBC8F",
  "Shared" = "#58B0E6"
)

flags_base <- phase_sig2 %>%
  dplyr::transmute(
    phase,
    scenario = "NAO_veg",
    sig_pure_A = !is.na(p_pure_A) & p_pure_A < 0.05,
    sig_pure_B = !is.na(p_pure_B) & p_pure_B < 0.05,
    sig_pure_C = NA,
    sig_full   = !is.na(p_full) & p_full < 0.05
  )
flags_tep <- phase_sig3 %>%
  dplyr::transmute(
    phase,
    scenario = "NAO_veg_tephra",
    sig_pure_A = !is.na(p_pure_A) & p_pure_A < 0.05,
    sig_pure_B = !is.na(p_pure_B) & p_pure_B < 0.05,
    sig_pure_C = !is.na(p_pure_C) & p_pure_C < 0.05,
    sig_full   = !is.na(p_full) & p_full < 0.05
  )
flags <- dplyr::bind_rows(flags_base, flags_tep)

plot_comps <- c(paste0("Pure ", A), paste0("Pure ", B), paste0("Pure ", C), "Shared")
vp_plot <- vp %>%
  dplyr::filter(component %in% plot_comps) %>%
  dplyr::filter(!is.na(value), value > 0) %>%
  dplyr::left_join(flags, by = c("phase", "scenario")) %>%
  dplyr::mutate(
    component_label = dplyr::recode(
      component,
      "Pure NAO_Median_Value" = "Pure NAO",
      "Pure estimate" = "Pure vegetation",
      "Pure tephra" = "Pure tephra"
    ),
    alpha_flag = dplyr::case_when(
      component_label == "Pure NAO" ~ ifelse(sig_pure_A, alpha_sig, alpha_ns),
      component_label == "Pure vegetation" ~ ifelse(sig_pure_B, alpha_sig, alpha_ns),
      component_label == "Pure tephra" ~ ifelse(sig_pure_C, alpha_sig, alpha_ns),
      component_label == "Shared" & shared_rule == "full_model" ~
        ifelse(sig_full, alpha_sig, alpha_ns),
      component_label == "Shared" ~ alpha_sig,
      TRUE ~ alpha_sig
    ),
    scenario = factor(
      scenario,
      levels = c("NAO_veg", "NAO_veg_tephra"),
      labels = c("NAO + vegetation", "NAO + vegetation + tephra")
    ),
    component_label = factor(
      component_label,
      levels = c("Pure NAO", "Pure vegetation", "Pure tephra", "Shared")
    ),
    phase = factor(phase, levels = paste("Phase", 1:5))
  )

# Percent labels (all positive fractions; NS fractions still labelled but faded)
label_df <- vp_plot %>%
  dplyr::filter(value > 0.02)

p <- ggplot2::ggplot(
  vp_plot,
  ggplot2::aes(x = phase, y = value, fill = component_label, alpha = alpha_flag)
) +
  ggplot2::geom_col(position = "stack", width = 0.7) +
  ggplot2::geom_text(
    data = label_df,
    ggplot2::aes(label = paste0(round(value * 100, 1), "%")),
    position = ggplot2::position_stack(vjust = 0.5),
    color = "black", size = 2.6, alpha = 1
  ) +
  ggplot2::facet_wrap(~scenario) +
  ggplot2::scale_alpha_identity() +
  ggplot2::scale_fill_manual(values = component_colors, name = NULL) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.caption = ggplot2::element_text(size = 8, hjust = 0, colour = "grey30")
  ) +
  ggplot2::labs(
    title = "A6: Lake-level varpart with lake-specific tephra (+ Condition(lake))",
    subtitle = "Faded fractions = non-significant (permutation p ≥ 0.05)",
    y = expression("Adjusted "*R^2),
    x = NULL,
    caption = paste0(
      "Tephra = 1 only for lake×30-yr bins within ±1 bin of that lake's ",
      "SUPPORTED/TENTATIVE ages (no regional tephra curve). ",
      "Condition(lake) partials out lake identity. n_perm = ", n_perm, "."
    )
  )

ggplot2::ggsave(
  "outputs/revision/figures/A6_tephra_varpart.png",
  p, width = 9.5, height = 5.5, dpi = 150, bg = "white"
)

message(
  "A6 complete (", nrow(tephra), " tephra ages; ",
  nrow(flagged), "/", nrow(lake_env), " lake×bins with tephra=1; ",
  "mean tephra prevalence = ", round(mean(lake_env$tephra), 3), ")"
)
