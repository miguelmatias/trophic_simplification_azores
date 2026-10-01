# A11 — Alternative Figure 5 with tephra overlay on all panels
# Rebuilds variance-partitioning panels (Fig. 5 a–c) and overlays the A6
# tephra signal: ±1×30-yr windows around lake-specific SUPPORTED/TENTATIVE
# ages (include_sensitivity). Overlay only — no data filtering / dropping.
source("scripts/revision/_bootstrap.R")
shared <- revision_bootstrap()

suppressPackageStartupMessages({
  library(patchwork)
  library(ggplot2)
  library(purrr)
  library(tidyr)
  library(dplyr)
})

message("=== A11 alternative Figure 5 with tephra ===")

# ---- Tephra windows (same logic as A06_volcanism_sensitivity.R) ----
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

bin_size <- 30
tephra_bins <- tephra %>%
  dplyr::mutate(
    age_bin = floor(age_ce / bin_size) * bin_size,
    influence_lo = age_bin - bin_size,
    influence_hi = age_bin + bin_size
  )

influence_ages <- sort(unique(c(
  tephra_bins$age_bin,
  tephra_bins$influence_lo,
  tephra_bins$influence_hi
)))
influence_ages <- influence_ages[is.finite(influence_ages)]

# Merge contiguous flagged bin centres into continuous shaded ribbons
# (each 30-yr bin centre spans [age − 15, age + 15])
merge_tephra_bands <- function(ages, half_width = bin_size / 2) {
  ages <- sort(unique(ages[is.finite(ages)]))
  if (length(ages) == 0) {
    return(tibble::tibble(xmin = numeric(), xmax = numeric()))
  }
  gaps <- diff(ages) > bin_size + 1e-8
  grp <- cumsum(c(TRUE, gaps))
  tibble::tibble(age = ages, grp = grp) %>%
    dplyr::group_by(grp) %>%
    dplyr::summarise(
      xmin = min(age) - half_width,
      xmax = max(age) + half_width,
      .groups = "drop"
    )
}

tephra_bands <- merge_tephra_bands(influence_ages) %>%
  dplyr::mutate(
    xmin = pmax(xmin, 0),
    xmax = pmin(xmax, 2010)
  ) %>%
  dplyr::filter(xmax > xmin, xmax >= 0, xmin <= 2010)

revision_write_csv(tephra_bands, "outputs/revision/A11_tephra_bands.csv")
revision_write_csv(
  tibble::tibble(age_ce = influence_ages[influence_ages >= 0 & influence_ages <= 2010]),
  "outputs/revision/A11_tephra_influence_ages.csv"
)

message(
  "Tephra ages used: ", nrow(tephra),
  "; influenced bin centres: ", length(influence_ages),
  "; merged bands: ", nrow(tephra_bands)
)

# ---- Rebuild Fig. 5 varpart data from shared joined table ----
joined_df_30yr <- shared$joined_df_30yr
stopifnot(!is.null(joined_df_30yr), nrow(joined_df_30yr) > 0)

selected_variables <- c("NAO_Median_Value", "estimate")
guild_cols <- c(
  "high_profile", "low_profile", "motile", "euplanktonic",
  "algivore", "detritivore", "plantivore", "predator"
)
A <- "NAO_Median_Value"
B <- "estimate"

# Faster than published 9999; overlay is the focus. Pattern should match closely.
n_perm <- 999
min_n <- 5
set.seed(1)

phase_list <- list(
  "Phase 1" = joined_df_30yr %>% dplyr::filter(age_ce < 750),
  "Phase 2" = joined_df_30yr %>% dplyr::filter(age_ce >= 750 & age_ce < 1050),
  "Phase 3" = joined_df_30yr %>% dplyr::filter(age_ce >= 1050 & age_ce < 1450),
  "Phase 4" = joined_df_30yr %>% dplyr::filter(age_ce >= 1450 & age_ce < 1750),
  "Phase 5" = joined_df_30yr %>% dplyr::filter(age_ce >= 1750)
)

phase_year_bounds <- tibble::tibble(
  phase = paste("Phase", 1:5),
  phase_lo = c(0, 750, 1050, 1450, 1750),
  phase_hi = c(750, 1050, 1450, 1750, 2010)
)

# Effects (adj. R²)
message("Computing phase varpart effects...")
r2_partitioned <- purrr::imap_dfr(phase_list, calculate_pure_shared_phase) %>%
  dplyr::mutate(
    component = factor(
      component,
      levels = c("Shared", "Pure Vegetation", "Pure Climate")
    ),
    phase_label = factor(
      phase,
      levels = paste("Phase", 1:5)
    )
  )

# Phase significance (fade NS fractions)
message("Computing phase significance (n_perm = ", n_perm, ")...")
phase_sig <- purrr::imap_dfr(phase_list, function(dat, phase_name) {
  sub <- dat %>%
    dplyr::select(age_ce, n_lakes, dplyr::all_of(guild_cols), dplyr::all_of(c(A, B))) %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(c(guild_cols, A, B)), ~ !is.na(.x)))

  n_bins <- nrow(sub)
  mean_n_lakes <- mean(sub$n_lakes, na.rm = TRUE)

  if (n_bins < min_n) {
    return(tibble::tibble(
      phase = phase_name,
      n_bins = n_bins,
      mean_n_lakes = mean_n_lakes,
      p_full = NA_real_,
      p_pure_A = NA_real_,
      p_pure_B = NA_real_
    ))
  }

  comm <- sub %>% dplyr::select(dplyr::all_of(guild_cols)) %>% as.matrix()
  env  <- sub %>% dplyr::select(dplyr::all_of(c(A, B))) %>% as.data.frame()

  mod_full <- vegan::rda(
    as.formula(paste("comm ~", paste(c(A, B), collapse = " + "))), data = env
  )
  mod_pure_A <- vegan::rda(
    as.formula(paste("comm ~", A, "+ Condition(", B, ")")), data = env
  )
  mod_pure_B <- vegan::rda(
    as.formula(paste("comm ~", B, "+ Condition(", A, ")")), data = env
  )

  tibble::tibble(
    phase = phase_name,
    n_bins = n_bins,
    mean_n_lakes = mean_n_lakes,
    p_full = safe_anova_p(mod_full, permutations = n_perm),
    p_pure_A = safe_anova_p(mod_pure_A, permutations = n_perm),
    p_pure_B = safe_anova_p(mod_pure_B, permutations = n_perm)
  )
})

alpha_sig <- 1.00
alpha_ns  <- 0.25

r2_partitioned_sig <- r2_partitioned %>%
  dplyr::left_join(
    phase_sig %>%
      dplyr::transmute(
        phase,
        sig_pure_A = p_pure_A < 0.05,
        sig_pure_B = p_pure_B < 0.05
      ),
    by = "phase"
  ) %>%
  dplyr::mutate(
    alpha_flag = dplyr::case_when(
      component == "Pure Climate"    ~ ifelse(sig_pure_A, alpha_sig, alpha_ns),
      component == "Pure Vegetation" ~ ifelse(sig_pure_B, alpha_sig, alpha_ns),
      component == "Shared"          ~ alpha_sig,
      TRUE ~ alpha_sig
    )
  )

# Moving-window effects
window_size <- 300
step_size   <- 30
min_year <- floor(min(joined_df_30yr$age_ce, na.rm = TRUE))
max_year <- ceiling(max(joined_df_30yr$age_ce, na.rm = TRUE))
window_starts <- seq(min_year, max_year - window_size, by = step_size)

message("Computing moving-window varpart effects (", length(window_starts), " windows)...")
r2_windowed <- purrr::map_dfr(
  window_starts,
  ~ calculate_pure_shared_window(joined_df_30yr, .x, window_size, selected_variables)
) %>%
  dplyr::mutate(
    midpoint = (window_start + window_end) / 2,
    component = factor(
      component,
      levels = c("Shared", "Pure Vegetation", "Pure Climate")
    )
  )

message("Computing moving-window significance (n_perm = ", n_perm, ")...")
window_sig <- purrr::map_dfr(
  window_starts,
  ~ calc_window_significance(
    joined_df_30yr,
    window_start = .x,
    window_size  = window_size,
    guild_cols   = guild_cols,
    A = A, B = B,
    min_n = min_n,
    n_perm = n_perm
  )
)

r2_windowed_sig <- r2_windowed %>%
  dplyr::left_join(
    window_sig %>%
      dplyr::transmute(
        window_start, window_end, midpoint,
        sig_pure_A = p_pure_A < 0.05,
        sig_pure_B = p_pure_B < 0.05
      ),
    by = c("window_start", "window_end", "midpoint")
  ) %>%
  dplyr::mutate(
    alpha_flag = dplyr::case_when(
      component == "Pure Climate"    ~ ifelse(sig_pure_A, alpha_sig, alpha_ns),
      component == "Pure Vegetation" ~ ifelse(sig_pure_B, alpha_sig, alpha_ns),
      component == "Shared"          ~ alpha_sig,
      TRUE ~ alpha_sig
    )
  )

# Panel (c) effect-size difference
effect_diff <- r2_windowed %>%
  dplyr::filter(component %in% c("Pure Vegetation", "Pure Climate")) %>%
  dplyr::select(window_start, window_end, midpoint, component, value) %>%
  tidyr::pivot_wider(names_from = component, values_from = value) %>%
  dplyr::left_join(
    window_sig %>% dplyr::select(window_start, window_end, midpoint, p_pure_A, p_pure_B),
    by = c("window_start", "window_end", "midpoint")
  ) %>%
  dplyr::mutate(
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
      direction == "veg"  & sig_veg  ~ 1.0,
      direction == "veg"  & !sig_veg ~ 0.25,
      direction == "clim" & sig_clim ~ 1.0,
      direction == "clim" & !sig_clim ~ 0.25,
      TRUE ~ 0.25
    )
  ) %>%
  dplyr::filter(!is.na(effect_group))

# Cache rebuild tables for traceability
revision_write_csv(r2_partitioned_sig, "outputs/revision/A11_fig5_phase_varpart.csv")
revision_write_csv(r2_windowed_sig, "outputs/revision/A11_fig5_window_varpart.csv")
revision_write_csv(effect_diff, "outputs/revision/A11_fig5_effect_diff.csv")

# ---- Plot styling (match Fig. 5) ----
component_colors <- c(
  "Pure Vegetation" = "#97D8C4",
  "Pure Climate" = "#F4B942",
  "Shared" = "#58B0E6"
)
effect_colors_split <- c(
  "VegChange > Climate (pre-750)"  = "#054A29",
  "VegChange > Climate (post-750)" = "#6DC6B6",
  "Climate > VegChange"            = "#F4B942"
)

tephra_fill <- "#8B7355"
tephra_alpha <- 0.18

x_min <- floor(min(joined_df_30yr$age_ce, na.rm = TRUE) / 100) * 100
x_max <- ceiling(max(joined_df_30yr$age_ce, na.rm = TRUE) / 100) * 100
common_breaks <- seq(x_min, x_max, by = 200)
common_minor  <- seq(x_min, x_max, by = 100)

base_size       <- 10
title_size      <- 11
axis_title_size <- 10
axis_text_size  <- 9
legend_title_sz <- 9
legend_text_sz  <- 8
fig_width  <- 4.5
fig_height <- 10.6

panel_theme <- ggplot2::theme_minimal(base_size = base_size) +
  ggplot2::theme(
    plot.title        = ggplot2::element_text(size = title_size, face = "bold", hjust = 0),
    axis.title        = ggplot2::element_text(size = axis_title_size),
    axis.text         = ggplot2::element_text(size = axis_text_size),
    legend.position   = "bottom",
    legend.title      = ggplot2::element_text(size = legend_title_sz),
    legend.text       = ggplot2::element_text(size = legend_text_sz),
    legend.key.height = grid::unit(3.0, "mm"),
    legend.key.width  = grid::unit(6.0, "mm"),
    plot.margin       = ggplot2::margin(3, 3, 3, 3)
  )

# Continuous-year tephra layer (panels b–c; also used on panel a timeline strip)
tephra_layer_year <- ggplot2::geom_rect(
  data = tephra_bands,
  ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
  inherit.aes = FALSE,
  fill = tephra_fill,
  alpha = tephra_alpha,
  colour = NA
)

# Panel (a): categorical phases — shade whole phases that overlap tephra windows
phase_tephra <- phase_year_bounds %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    has_tephra = any(tephra_bands$xmin < phase_hi & tephra_bands$xmax > phase_lo)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(has_tephra) %>%
  dplyr::mutate(
    phase_label = factor(phase, levels = paste("Phase", 1:5)),
    xmin = as.numeric(phase_label) - 0.45,
    xmax = as.numeric(phase_label) + 0.45
  )

tephra_layer_phase <- ggplot2::geom_rect(
  data = phase_tephra,
  ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
  inherit.aes = FALSE,
  fill = tephra_fill,
  alpha = tephra_alpha,
  colour = NA
)

tephra_legend_df <- data.frame(
  x = NA_real_, y = NA_real_,
  band = "Tephra (±1×30-yr)"
)

# ---- Panel (a): Historical phases ----
p_a <- ggplot2::ggplot(
  r2_partitioned_sig,
  ggplot2::aes(x = phase_label, y = value, fill = component, alpha = alpha_flag)
) +
  tephra_layer_phase +
  ggplot2::geom_col(width = 0.6) +
  ggplot2::geom_text(
    data = r2_partitioned_sig %>% dplyr::filter(value > 0),
    ggplot2::aes(label = paste0(round(value * 100, 1), "%")),
    position = ggplot2::position_stack(vjust = 0.5),
    color = "black", size = 2.4
  ) +
  # Invisible point for tephra legend key (separate fill scale via new_scale if needed;
  # here we append tephra as a fill level via a dummy geom + manual values)
  ggplot2::geom_point(
    data = tephra_legend_df,
    ggplot2::aes(x = x, y = y, fill = band),
    shape = 22, size = 4, colour = NA, alpha = 0, inherit.aes = FALSE
  ) +
  ggplot2::scale_alpha_identity() +
  ggplot2::scale_fill_manual(
    name = "Variance Component",
    values = c(component_colors, "Tephra (±1×30-yr)" = tephra_fill),
    breaks = c("Pure Climate", "Pure Vegetation", "Shared", "Tephra (±1×30-yr)"),
    guide = ggplot2::guide_legend(
      nrow = 2, byrow = TRUE,
      override.aes = list(
        alpha = c(1, 1, 1, 0.55),
        shape = c(NA, NA, NA, 22),
        size = c(NA, NA, NA, 5)
      )
    )
  ) +
  ggplot2::scale_x_discrete(labels = paste("Phase", 1:5), drop = FALSE) +
  ggplot2::labs(
    x = NULL,
    y = expression("Adjusted "*R^2*" (Variance explained)"),
    title = "(a) Historical phases"
  ) +
  panel_theme

# ---- Panel (b): Moving window ----
p_b <- ggplot2::ggplot() +
  tephra_layer_year +
  ggplot2::geom_col(
    data = r2_windowed_sig,
    ggplot2::aes(x = midpoint, y = value, fill = component, alpha = alpha_flag),
    width = step_size
  ) +
  ggplot2::geom_point(
    data = tephra_legend_df,
    ggplot2::aes(x = x, y = y, fill = band),
    shape = 22, size = 4, colour = NA, alpha = 0, inherit.aes = FALSE
  ) +
  ggplot2::scale_alpha_identity() +
  ggplot2::scale_fill_manual(
    name = "Variance Component",
    values = c(component_colors, "Tephra (±1×30-yr)" = tephra_fill),
    breaks = c("Pure Climate", "Pure Vegetation", "Shared", "Tephra (±1×30-yr)"),
    guide = ggplot2::guide_legend(
      nrow = 2, byrow = TRUE,
      override.aes = list(
        alpha = c(1, 1, 1, 0.55),
        shape = c(NA, NA, NA, 22),
        size = c(NA, NA, NA, 5)
      )
    )
  ) +
  ggplot2::scale_x_continuous(
    limits = c(x_min, x_max),
    breaks = common_breaks,
    minor_breaks = common_minor,
    expand = c(0, 0)
  ) +
  ggplot2::labs(
    x = "Year (CE)",
    y = expression("Adjusted "*R^2*" (Variance explained)"),
    title = "(b) Moving window"
  ) +
  panel_theme

# ---- Panel (c): Climate vs VegChange ----
p_c <- ggplot2::ggplot(
  effect_diff,
  ggplot2::aes(x = midpoint, y = effect_diff, fill = effect_group, alpha = alpha_flag)
) +
  tephra_layer_year +
  ggplot2::geom_col(width = step_size) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  ggplot2::geom_point(
    data = tephra_legend_df,
    ggplot2::aes(x = x, y = y, fill = band),
    shape = 22, size = 4, colour = NA, alpha = 0, inherit.aes = FALSE
  ) +
  ggplot2::scale_alpha_identity() +
  ggplot2::scale_fill_manual(
    name = "Dominant Driver",
    values = c(effect_colors_split, "Tephra (±1×30-yr)" = tephra_fill),
    breaks = c(
      "Climate > VegChange",
      "VegChange > Climate (post-750)",
      "VegChange > Climate (pre-750)",
      "Tephra (±1×30-yr)"
    ),
    drop = FALSE,
    guide = ggplot2::guide_legend(
      nrow = 2, byrow = TRUE,
      override.aes = list(
        alpha = c(1, 1, 1, 0.55),
        shape = c(NA, NA, NA, 22),
        size = c(NA, NA, NA, 5)
      )
    )
  ) +
  ggplot2::scale_x_continuous(
    limits = c(x_min, x_max),
    breaks = common_breaks,
    minor_breaks = common_minor,
    expand = c(0, 0)
  ) +
  ggplot2::labs(
    x = "Year (CE)",
    y = "Effect size (Vegetation − Climate)",
    title = "(c) Climate vs VegChange",
    caption = paste0(
      "Bands = lake-specific tephra ages (±1×30-yr bin around SUPPORTED/TENTATIVE ",
      "events; same coding as A6). Panel (a): phases overlapping any tephra window ",
      "are shaded; (b–c): vertical year bands. Overlay only — no bins dropped. ",
      "Significance fading uses n_perm = ", n_perm, " (published Fig. 5 used 9999)."
    )
  ) +
  panel_theme +
  ggplot2::theme(
    plot.caption = ggplot2::element_text(
      size = 7, hjust = 0, colour = "grey30", lineheight = 1.1
    )
  )

alt_fig5 <- (p_a / p_b / p_c) +
  patchwork::plot_layout(guides = "keep", heights = c(1, 1.2, 1.15))

out_png <- "outputs/revision/figures/A11_alt_figure5_with_tephra.png"
out_pdf <- "outputs/revision/figures/A11_alt_figure5_with_tephra.pdf"

ggplot2::ggsave(
  out_png, alt_fig5,
  width = fig_width, height = fig_height, units = "in", dpi = 300, bg = "white"
)

ggplot2::ggsave(
  out_pdf, alt_fig5,
  width = fig_width, height = fig_height, units = "in", device = "pdf", bg = "white"
)

if (!file.exists(out_pdf)) {
  warning("PDF was not written: ", out_pdf)
} else {
  message("Saved ", out_png, " and ", out_pdf)
}
message("A11 complete")
