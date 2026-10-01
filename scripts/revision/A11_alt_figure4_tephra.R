# A11 — Alternative Figure 4 with tephra overlay on all panels
# Rebuilds regional standardised guild richness (Fig. 4 a–b) and adds an
# overall producers-vs-consumers panel (c). Tephra windows follow A6:
# ±1×30-yr bin around lake-specific SUPPORTED/TENTATIVE ages
# (include_sensitivity). Overlay only — no data filtering.
source("scripts/revision/_bootstrap.R")
shared <- revision_bootstrap()

suppressPackageStartupMessages({
  library(patchwork)
  library(scico)
  library(viridis)
  library(ggplot2)
})

message("=== A11 alternative Figure 4 with tephra ===")

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

# ---- Fig. 4 data ----
prod_df <- shared$prod_glob_div_plt_data
cons_df <- shared$cons_glob_div_plt_data
if (is.null(prod_df) || is.null(cons_df)) {
  e <- new.env(parent = emptyenv())
  load("data/clean_source_data_files.RData", envir = e)
  prod_df <- get("prod_glob_div_plt_data", envir = e)
  cons_df <- get("cons_glob_div_plt_data", envir = e)
}

cons_labels <- c(
  predator = "Predators",
  detritivore = "Detritivores",
  algivore = "Algivores",
  plantivore = "Plantivores"
)
prod_labels <- c(
  high_profile = "High-profile",
  low_profile = "Low-profile",
  motile = "Motile",
  euplanktonic = "Euplanktonic"
)

axis_title_size <- 13
axis_text_size <- 11
tephra_fill <- "#8B7355"  # warm ash / taupe
tephra_alpha <- 0.18

# Fixed fill (not mapped) so geom_smooth SE ribbons keep their own fill
tephra_layer <- ggplot2::geom_rect(
  data = tephra_bands,
  ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
  inherit.aes = FALSE,
  fill = tephra_fill,
  alpha = tephra_alpha,
  colour = NA
)

# Dummy data for a single shared tephra legend key
tephra_legend_df <- data.frame(
  x = NA_real_, y = NA_real_,
  band = "Tephra (±1×30-yr)"
)

panel_theme <- ggplot2::theme_minimal() +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(size = axis_title_size),
    axis.text = ggplot2::element_text(size = axis_text_size),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.direction = "horizontal",
    legend.title = ggplot2::element_text(size = 11),
    legend.text = ggplot2::element_text(size = 10),
    legend.key.width = grid::unit(1.2, "lines"),
    legend.key.height = grid::unit(0.85, "lines"),
    plot.margin = ggplot2::margin(5.5, 5.5, 8, 5.5)
  )

# Panel (a): producers by guild — matches published Fig. 4a
p_a <- ggplot2::ggplot(
  prod_df,
  ggplot2::aes(age_ce, standardized_species_richness)
) +
  tephra_layer +
  ggplot2::geom_smooth(ggplot2::aes(color = fgroup), se = TRUE) +
  ggplot2::geom_smooth(
    color = "red", linetype = "dotted", alpha = 0.7, se = FALSE
  ) +
  # Invisible point for tephra legend entry (fill aesthetic isolated)
  ggplot2::geom_point(
    data = tephra_legend_df,
    ggplot2::aes(x = x, y = y, fill = band),
    shape = 22, size = 4, colour = NA, alpha = 0, inherit.aes = FALSE
  ) +
  ggplot2::scale_color_viridis_d(
    name = "Producers",
    option = "D",
    direction = 1,
    alpha = 0.5,
    labels = prod_labels
  ) +
  ggplot2::scale_fill_manual(
    name = NULL,
    values = c("Tephra (±1×30-yr)" = tephra_fill),
    guide = ggplot2::guide_legend(
      override.aes = list(alpha = 0.55, size = 5, shape = 22, colour = NA)
    )
  ) +
  ggplot2::labs(
    x = "Age (CE)",
    y = "Standardized Species Richness (0–1)",
    tag = "a"
  ) +
  panel_theme +
  ggplot2::guides(
    color = ggplot2::guide_legend(nrow = 1, byrow = TRUE, order = 1)
  ) +
  ggplot2::lims(x = c(0, 2010))

# Panel (b): consumers by guild — matches published Fig. 4b
p_b <- ggplot2::ggplot(
  cons_df,
  ggplot2::aes(age_ce, standardized_species_richness)
) +
  tephra_layer +
  ggplot2::geom_smooth(ggplot2::aes(color = fgroup), se = TRUE) +
  ggplot2::geom_smooth(
    color = "red", linetype = "dotted", alpha = 0.7, se = FALSE
  ) +
  ggplot2::geom_point(
    data = tephra_legend_df,
    ggplot2::aes(x = x, y = y, fill = band),
    shape = 22, size = 4, colour = NA, alpha = 0, inherit.aes = FALSE
  ) +
  scico::scale_color_scico_d(
    name = "Consumers",
    palette = "lajolla",
    direction = -1,
    labels = cons_labels
  ) +
  ggplot2::scale_fill_manual(
    name = NULL,
    values = c("Tephra (±1×30-yr)" = tephra_fill),
    guide = ggplot2::guide_legend(
      override.aes = list(alpha = 0.55, size = 5, shape = 22, colour = NA)
    )
  ) +
  ggplot2::labs(
    x = "Age (CE)",
    y = "Standardized Species Richness (0–1)",
    tag = "b"
  ) +
  panel_theme +
  ggplot2::guides(
    color = ggplot2::guide_legend(nrow = 1, byrow = TRUE, order = 1)
  ) +
  ggplot2::lims(x = c(0, 2010))

# Panel (c): overall producers vs consumers (aggregate trajectories)
# Published Fig. 4 has only a–b; this panel shows the overall smooths
# as two series so the alt figure has three panels (user request a–c).
overall_df <- dplyr::bind_rows(
  prod_df %>% dplyr::mutate(trophic = "Producers"),
  cons_df %>% dplyr::mutate(trophic = "Consumers")
) %>%
  dplyr::mutate(
    trophic = factor(trophic, levels = c("Producers", "Consumers"))
  )

p_c <- ggplot2::ggplot(
  overall_df,
  ggplot2::aes(age_ce, standardized_species_richness)
) +
  tephra_layer +
  ggplot2::geom_smooth(
    ggplot2::aes(color = trophic),
    se = TRUE,
    linewidth = 1
  ) +
  ggplot2::geom_point(
    data = tephra_legend_df,
    ggplot2::aes(x = x, y = y, fill = band),
    shape = 22, size = 4, colour = NA, alpha = 0, inherit.aes = FALSE
  ) +
  ggplot2::scale_color_manual(
    name = "Overall",
    values = c(Producers = "#FDE725", Consumers = "#440154")
  ) +
  ggplot2::scale_fill_manual(
    name = NULL,
    values = c("Tephra (±1×30-yr)" = tephra_fill),
    guide = ggplot2::guide_legend(
      override.aes = list(alpha = 0.55, size = 5, shape = 22, colour = NA)
    )
  ) +
  ggplot2::labs(
    x = "Age (CE)",
    y = "Standardized Species Richness (0–1)",
    tag = "c",
    caption = paste0(
      "Bands = lake-specific tephra ages (±1×30-yr bin around SUPPORTED/TENTATIVE ",
      "events; same coding as A6). Panels a–b match published Fig. 4; ",
      "panel c = overall producer/consumer smooths (not in published Fig. 4)."
    )
  ) +
  panel_theme +
  ggplot2::theme(
    plot.caption = ggplot2::element_text(
      size = 8.5, hjust = 0, colour = "grey30"
    )
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(nrow = 1, byrow = TRUE, order = 1)
  ) +
  ggplot2::lims(x = c(0, 2010))

alt_fig4 <- (p_a / p_b / p_c) +
  patchwork::plot_layout(heights = c(1, 1, 1))

out_png <- "outputs/revision/figures/A11_alt_figure4_with_tephra.png"
out_pdf <- "outputs/revision/figures/A11_alt_figure4_with_tephra.pdf"

ggplot2::ggsave(
  out_png, alt_fig4,
  width = 7.2, height = 14.5, units = "in", dpi = 300
)

# PDF via default device (cairo_pdf can fail in some headless/sandbox runs)
ggplot2::ggsave(
  out_pdf, alt_fig4,
  width = 7.2, height = 14.5, units = "in", device = "pdf"
)

if (!file.exists(out_pdf)) {
  warning("PDF was not written: ", out_pdf)
} else {
  message("Saved ", out_png, " and ", out_pdf)
}
message("A11 complete")
