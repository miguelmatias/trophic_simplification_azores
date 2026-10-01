# Revision-analysis helpers (A1–A9)
# Sourced by scripts under scripts/revision/

revision_guild_cols <- c(
  "high_profile", "low_profile", "motile", "euplanktonic",
  "algivore", "detritivore", "plantivore", "predator"
)

revision_phase_breaks <- list(
  "Phase 1" = c(-Inf, 750),
  "Phase 2" = c(750, 1050),
  "Phase 3" = c(1050, 1450),
  "Phase 4" = c(1450, 1750),
  "Phase 5" = c(1750, Inf)
)

revision_ensure_dirs <- function(root = ".") {
  dirs <- c(
    file.path(root, "outputs/revision"),
    file.path(root, "outputs/revision/figures"),
    file.path(root, "data/revision")
  )
  lapply(dirs, dir.create, showWarnings = FALSE, recursive = TRUE)
  invisible(TRUE)
}

revision_meta_cols <- function() {
  c("lake", "core_depth_id", "depth_top", "depth_bottom", "dec_depth", "age_ce", "age_bp")
}

revision_taxon_cols <- function(df) {
  meta <- revision_meta_cols()
  cols <- setdiff(names(df), meta)
  cols[vapply(df[cols], is.numeric, logical(1))]
}

#' Wide count matrix + ages for one lake list
revision_bind_wide_with_ages <- function(wide_list, codes_list) {
  purrr::map2_dfr(wide_list, codes_list, function(wide, codes) {
    tax <- revision_taxon_cols(wide)
    ages <- codes %>%
      dplyr::select(dplyr::any_of(c("lake", "core_depth_id", "age_ce", "dec_depth")))
    wide %>%
      dplyr::select(dplyr::any_of(c("lake", "core_depth_id")), dplyr::all_of(tax)) %>%
      dplyr::left_join(ages, by = c("lake", "core_depth_id")) %>%
      dplyr::mutate(
        total_count = rowSums(dplyr::across(dplyr::all_of(tax)), na.rm = TRUE),
        age_ce = as.numeric(age_ce)
      )
  })
}

#' Pool counts within lake × 30-yr bins
revision_pool_bins <- function(wide_ages, bin_size = 30) {
  tax <- revision_taxon_cols(wide_ages)
  wide_ages %>%
    dplyr::filter(is.finite(age_ce), age_ce > 0) %>%
    dplyr::mutate(age_bin = floor(age_ce / bin_size) * bin_size) %>%
    dplyr::group_by(lake, age_bin) %>%
    dplyr::summarise(
      n_samples = dplyr::n(),
      total_count = sum(total_count, na.rm = TRUE),
      dplyr::across(dplyr::all_of(tax), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::rename(age_ce = age_bin)
}

#' Choose rarefaction depth retaining >= target fraction of bins
revision_choose_rarefy_n <- function(totals, target_keep = 0.70, candidates = NULL) {
  totals <- totals[is.finite(totals) & totals > 0]
  if (length(totals) == 0) return(NA_integer_)
  if (is.null(candidates)) {
    candidates <- unique(c(10L, 20L, 30L, 40L, 50L, 100L, 200L, 300L, 400L, 500L))
    candidates <- candidates[candidates <= max(totals)]
  }
  keep <- vapply(candidates, function(n) mean(totals >= n), numeric(1))
  ok <- candidates[keep >= target_keep]
  if (length(ok) == 0) {
    ok <- candidates[keep >= 0.5]
  }
  if (length(ok) == 0) return(as.integer(floor(stats::median(totals))))
  as.integer(max(ok))
}

#' Rarefy one count matrix row to n individuals; return richness
revision_rarefy_richness <- function(counts, n, n_rep = 20L, seed = 1L) {
  counts <- as.numeric(counts)
  counts[!is.finite(counts)] <- 0
  counts[counts < 0] <- 0
  counts <- as.integer(round(counts))
  total <- sum(counts)
  if (!is.finite(n) || n < 1 || total < n) return(NA_real_)
  n <- as.integer(n)
  if (requireNamespace("vegan", quietly = TRUE)) {
    set.seed(seed)
    vals <- replicate(n_rep, {
      rare <- vegan::rrarefy(matrix(counts, nrow = 1), sample = n)
      sum(rare > 0)
    })
    return(mean(vals))
  }
  set.seed(seed)
  pool <- rep.int(seq_along(counts), times = counts)
  if (length(pool) < n) return(NA_real_)
  mean(replicate(n_rep, length(unique(sample(pool, n)))))
}

revision_extract_genus <- function(taxon_names) {
  gsub("^([A-Za-z]+).*", "\\1", taxon_names)
}

#' Aggregate diatom wide matrices to genus level
revision_aggregate_diatoms_to_genus <- function(wide_list, guild_lookup = NULL) {
  if (is.null(guild_lookup)) {
    if (requireNamespace("readr", quietly = TRUE)) {
      guild_lookup <- as.data.frame(
        readr::read_csv("data/table_guilds_diatoms.csv", show_col_types = FALSE),
        stringsAsFactors = FALSE
      )
    } else {
      raw <- readBin("data/table_guilds_diatoms.csv", what = "raw", n = file.info("data/table_guilds_diatoms.csv")$size)
      if (length(raw) >= 3 && identical(raw[1:3], as.raw(c(0xef, 0xbb, 0xbf)))) {
        raw <- raw[-(1:3)]
      }
      txt <- rawToChar(raw)
      Encoding(txt) <- "UTF-8"
      guild_lookup <- utils::read.csv(text = txt, stringsAsFactors = FALSE, check.names = FALSE)
    }
  }
  names(guild_lookup)[1] <- "full_name"
  if (!"genus" %in% names(guild_lookup)) {
    stop("guild lookup must contain a 'genus' column")
  }
  lookup <- guild_lookup %>%
    dplyr::mutate(
      full_name = as.character(full_name),
      genus = as.character(genus)
    )

  purrr::imap(wide_list, function(df, lake_name) {
    tax <- revision_taxon_cols(df)
    meta <- df %>% dplyr::select(dplyr::any_of(c("lake", "core_depth_id")))
    mat <- as.matrix(df[tax])
    mat[!is.finite(mat)] <- 0

    map_df <- tibble::tibble(
      taxon = tax,
      genus_guess = revision_extract_genus(tax)
    ) %>%
      dplyr::left_join(
        lookup %>% dplyr::select(full_name, genus),
        by = c("taxon" = "full_name")
      ) %>%
      dplyr::mutate(genus = dplyr::coalesce(genus, genus_guess))

    genera <- unique(map_df$genus)
    out <- sapply(genera, function(g) {
      cols <- map_df$taxon[map_df$genus == g]
      rowSums(mat[, cols, drop = FALSE], na.rm = TRUE)
    })
    if (is.null(dim(out))) {
      out <- matrix(out, ncol = 1, dimnames = list(NULL, genera))
    }
    dplyr::bind_cols(meta, as.data.frame(out, check.names = FALSE))
  })
}

#' DCA1 scores for a list of wide abundance tables
revision_compute_dca1_list <- function(wide_list, codes_list, flip_lakes = character()) {
  names_w <- names(wide_list)
  if (is.null(names_w)) names_w <- names(codes_list)

  purrr::map2_dfr(wide_list, codes_list, function(wide, codes) {
    tax <- revision_taxon_cols(wide)
    mat <- wide %>%
      dplyr::select(dplyr::all_of(tax)) %>%
      replace(is.na(.), 0) %>%
      dplyr::select(where(~ sum(.x) > 0))
    if (ncol(mat) < 2 || nrow(mat) < 3) {
      return(tibble::tibble(
        lake = codes$lake[seq_len(nrow(codes))],
        core_depth_id = codes$core_depth_id,
        age_ce = as.numeric(codes$age_ce),
        DCA1 = NA_real_
      ))
    }
    hell <- vegan::decostand(mat, method = "hellinger")
    dca <- vegan::decorana(hell)
    sc <- as.data.frame(vegan::scores(dca, display = "sites", choices = 1))
    names(sc) <- "DCA1"
    out <- dplyr::bind_cols(
      codes %>% dplyr::select(dplyr::any_of(c("lake", "core_depth_id", "age_ce"))),
      sc
    ) %>%
      dplyr::mutate(age_ce = as.numeric(age_ce), DCA1 = as.numeric(DCA1))
    lk <- unique(out$lake)[1]
    if (lk %in% flip_lakes) out$DCA1 <- -out$DCA1
    out
  })
}

#' Build regional vegetation DCA1 HGAM smooth (yearly)
revision_build_vegetation_series <- function(pollen_df) {
  pollen_df <- pollen_df %>%
    dplyr::rename(age_bp = age) %>%
    dplyr::mutate(
      age_ce = 1950 - age_bp,
      ID = paste0(sitename, "_", age_ce)
    )

  exclude_sites <- c("Lomba", "Rasa")
  df_wk_dca <- pollen_df %>%
    dplyr::filter(
      !ecologicalgroup %in% c("UNID"),
      age_ce > 0,
      !sitename %in% exclude_sites
    )

  results_list <- lapply(sort(unique(df_wk_dca$sitename)), function(site) {
    df.sp <- df_wk_dca %>%
      dplyr::filter(sitename == site) %>%
      dplyr::group_by(ID, taxa) %>%
      dplyr::summarise(value = sum(value), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = taxa, values_from = value, values_fill = 0) %>%
      tibble::column_to_rownames("ID") %>%
      dplyr::select(where(~ any(. != 0)))
    if (nrow(df.sp) < 3 || ncol(df.sp) < 2) return(NULL)
    df.spe <- vegan::decostand(df.sp, method = "hellinger")
    vegan::decorana(df.spe) %>%
      vegan::scores(choices = 1) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("ID")
  })
  results_list <- Filter(Negate(is.null), results_list)
  DCA.res <- dplyr::bind_rows(results_list)

  score_cols <- setdiff(names(DCA.res), "ID")
  if (!"DCA1" %in% names(DCA.res)) {
    names(DCA.res)[names(DCA.res) == score_cols[1]] <- "DCA1"
  }

  DCA.res <- DCA.res %>%
    tidyr::separate(ID, into = c("sitename", "age_ce"), sep = "_", convert = TRUE) %>%
    dplyr::arrange(sitename, age_ce) %>%
    dplyr::mutate(
      DCA1 = dplyr::case_when(
        sitename == "Azul" ~ -DCA1,
        sitename == "Gingal" ~ -DCA1,
        TRUE ~ DCA1
      ),
      sitename = as.factor(sitename)
    )

  mod_to <- mgcv::gam(
    DCA1 ~ s(age_ce, bs = "tp", k = 25) +
      s(age_ce, by = sitename, bs = "tp", k = 15, m = 1) +
      s(sitename, bs = "re"),
    data = DCA.res,
    method = "REML"
  )

  interp_years <- 0:2010
  ref_site <- levels(DCA.res$sitename)[1]
  newdata_veg <- tibble::tibble(
    age_ce = interp_years,
    sitename = factor(ref_site, levels = levels(DCA.res$sitename))
  )
  pr_terms <- predict(mod_to, newdata = newdata_veg, type = "terms", se.fit = TRUE)
  term_names <- colnames(pr_terms$fit)
  idx <- which(term_names == "s(age_ce)")
  if (length(idx) == 0) idx <- grep("^s\\(age_ce\\)$", term_names)
  est <- as.numeric(pr_terms$fit[, idx])
  se <- as.numeric(pr_terms$se.fit[, idx])
  zcrit <- qnorm(0.975)

  veg_yearly <- tibble::tibble(
    age_ce = as.integer(interp_years),
    estimate = est,
    lower_ci = est - zcrit * se,
    upper_ci = est + zcrit * se
  )

  list(dca = DCA.res, mod = mod_to, veg_yearly = veg_yearly)
}

#' Arboreal percentage by site × age
revision_arboreal_series <- function(pollen_df) {
  pollen_df %>%
    dplyr::rename(age_bp = age) %>%
    dplyr::mutate(age_ce = 1950 - age_bp) %>%
    dplyr::filter(age_ce > 0, !ecologicalgroup %in% c("UNID", "AQVP", "AQBR")) %>%
    dplyr::group_by(sitename, age_ce, group) %>%
    dplyr::summarise(value = sum(value), .groups = "drop") %>%
    dplyr::group_by(sitename, age_ce) %>%
    dplyr::mutate(pct = 100 * value / sum(value)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(group == "Tree") %>%
    dplyr::select(sitename, age_ce, arboreal_pct = pct)
}

#' First presence of indicator taxa
revision_indicator_onsets <- function(pollen_df,
                                      taxa_cereals = c("Cerealia", "Triticum", "Secale", "Zea mays"),
                                      taxa_plantago = c(
                                        "Plantago", "Plantago lanceolata",
                                        "Plantago major", "Plantago coronopus"
                                      )) {
  pollen_df <- pollen_df %>%
    dplyr::rename(age_bp = age) %>%
    dplyr::mutate(age_ce = 1950 - age_bp)

  one_cat <- function(taxa, label) {
    pollen_df %>%
      dplyr::filter(taxa %in% .env$taxa, value > 0, age_ce > 0) %>%
      dplyr::group_by(sitename, age_ce) %>%
      dplyr::summarise(count = sum(value), .groups = "drop") %>%
      dplyr::arrange(sitename, age_ce) %>%
      dplyr::group_by(sitename) %>%
      dplyr::summarise(
        category = label,
        first_presence_ce = dplyr::first(age_ce),
        n_positive_samples = dplyr::n(),
        .groups = "drop"
      )
  }

  dplyr::bind_rows(
    one_cat(taxa_cereals, "Cereals"),
    one_cat(taxa_plantago, "Plantago")
  )
}

#' Build 30-yr joined biology + NAO + vegetation table
revision_build_joined_30yr <- function(norm_abund_guilds,
                                       nao_df,
                                       veg_yearly,
                                       bin_size = 30) {
  guild_cols <- revision_guild_cols
  bio_lake_30 <- norm_abund_guilds %>%
    dplyr::select(lake, age_ce, dplyr::all_of(guild_cols)) %>%
    dplyr::filter(age_ce > 0) %>%
    dplyr::mutate(
      age_ce = as.numeric(age_ce),
      age_bin = floor(age_ce / bin_size) * bin_size
    ) %>%
    dplyr::group_by(lake, age_bin) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(guild_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::rename(age_ce = age_bin)

  bio_reg_30 <- bio_lake_30 %>%
    dplyr::group_by(age_ce) %>%
    dplyr::summarise(
      n_lakes = dplyr::n_distinct(lake),
      dplyr::across(dplyr::all_of(guild_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::arrange(age_ce)

  nao <- nao_df
  if (!"NAO_Median_Value" %in% names(nao)) {
    names(nao)[2] <- "age_ce"
    names(nao)[5] <- "NAO_Median_Value"
  }
  nao <- nao %>%
    dplyr::mutate(
      age_ce = as.numeric(age_ce),
      NAO_Median_Value = as.numeric(NAO_Median_Value)
    ) %>%
    dplyr::filter(age_ce >= 0) %>%
    dplyr::select(age_ce, NAO_Median_Value)

  # local bin30 that does not rely on env bin_size
  bin_local <- function(df, cols) {
    df %>%
      dplyr::filter(!is.na(age_ce)) %>%
      dplyr::mutate(age_bin = floor(age_ce / bin_size) * bin_size) %>%
      dplyr::group_by(age_bin) %>%
      dplyr::summarise(
        dplyr::across(dplyr::all_of(cols), ~ mean(.x, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      dplyr::rename(age_ce = age_bin) %>%
      dplyr::arrange(age_ce)
  }

  nao_30 <- bin_local(nao, "NAO_Median_Value")
  veg_30 <- bin_local(veg_yearly, c("estimate", "lower_ci", "upper_ci"))

  joined <- bio_reg_30 %>%
    dplyr::left_join(nao_30, by = "age_ce") %>%
    dplyr::left_join(veg_30 %>% dplyr::select(age_ce, estimate), by = "age_ce")

  list(
    bio_lake_30 = bio_lake_30,
    bio_reg_30 = bio_reg_30,
    nao_30 = nao_30,
    veg_30 = veg_30,
    joined_df_30yr = joined
  )
}

#' Island-balanced regional mean of guild abundances
revision_island_balanced_bio_reg <- function(bio_lake_30, lake_island) {
  guild_cols <- revision_guild_cols
  bio_lake_30 %>%
    dplyr::left_join(lake_island, by = "lake") %>%
    dplyr::filter(!is.na(island)) %>%
    dplyr::group_by(age_ce, island) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(guild_cols), ~ mean(.x, na.rm = TRUE)),
      n_lakes_island = dplyr::n_distinct(lake),
      .groups = "drop"
    ) %>%
    dplyr::group_by(age_ce) %>%
    dplyr::summarise(
      n_islands = dplyr::n_distinct(island),
      n_lakes = sum(n_lakes_island),
      dplyr::across(dplyr::all_of(guild_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )
}

revision_cor_pair <- function(x, y, method = c("pearson", "spearman")) {
  method <- match.arg(method)
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4) {
    return(tibble::tibble(method = method, n = sum(ok), estimate = NA_real_, p.value = NA_real_))
  }
  ct <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = method))
  tibble::tibble(
    method = method,
    n = sum(ok),
    estimate = unname(ct$estimate),
    p.value = ct$p.value
  )
}

revision_normalize_lake_names <- function(x) {
  x <- as.character(x)
  x <- iconv(x, from = "", to = "UTF-8", sub = "")
  x <- trimws(x)
  # ASCII-safe replacements for known accents in this project
  x <- gsub("Caldeir.?o", "Caldeirao", x, perl = TRUE)
  x <- gsub("Emp\\.\\s*Norte", "Empadadas Norte", x)
  x <- gsub("^Emp Norte$", "Empadadas Norte", x)
  x
}

revision_load_lake_meta <- function(path = "data/table_lake_metadata.csv") {
  if (requireNamespace("readr", quietly = TRUE)) {
    df <- readr::read_csv(path, show_col_types = FALSE, locale = readr::locale(encoding = "UTF-8"))
    df <- as.data.frame(df, stringsAsFactors = FALSE)
  } else {
    raw <- readBin(path, what = "raw", n = file.info(path)$size)
    # strip UTF-8 BOM if present
    if (length(raw) >= 3 && raw[1] == as.raw(0xef) && raw[2] == as.raw(0xbb) && raw[3] == as.raw(0xbf)) {
      raw <- raw[-(1:3)]
    }
    txt <- rawToChar(raw)
    Encoding(txt) <- "UTF-8"
    df <- utils::read.csv(text = txt, stringsAsFactors = FALSE, check.names = FALSE)
  }
  names(df) <- gsub("^\ufeff", "", names(df))
  names(df) <- trimws(names(df))
  df %>%
    dplyr::mutate(
      lake = revision_normalize_lake_names(lake),
      island = trimws(as.character(island)),
      alt = as.numeric(alt),
      area = as.numeric(area),
      zmax = as.numeric(zmax)
    )
}

#' Viridis CTS palette matching main_script (CTS1 = yellow)
revision_cts_colours <- function(k = 5L) {
  cols <- viridis::viridis(as.integer(k), direction = -1)
  stats::setNames(cols, paste0("CTS", seq_len(k)))
}

#' Assign AMD clusters at fixed k without lumping
#' Clusters are relabelled by mean euplanctonic relative abundance
#' (highest → CTS1), matching the original MS yellow/euplanctonic CTS1 logic.
revision_assign_cts <- function(norm_abund, k = 6L, iterations = 400L) {
  mat <- norm_abund %>%
    dplyr::select(dplyr::all_of(revision_guild_cols)) %>%
    replace(is.na(.), 0)
  cl <- getAMDclusters(mat, .iterations = iterations, .opt_num_clusts = k)
  assign_vec <- as.integer(cl[[1]])
  out <- norm_abund %>%
    dplyr::mutate(amd_raw = factor(assign_vec))

  key <- out %>%
    dplyr::group_by(amd_raw) %>%
    dplyr::summarise(
      mean_euplanktonic = mean(euplanktonic, na.rm = TRUE),
      mean_consumer = mean(algivore + detritivore + plantivore + predator, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    # Highest euplanctonic → CTS1 (yellow); consumers break remaining ties
    dplyr::arrange(dplyr::desc(mean_euplanktonic), mean_consumer) %>%
    dplyr::mutate(rank = dplyr::row_number())

  map_rank <- stats::setNames(key$rank, as.character(key$amd_raw))
  out %>%
    dplyr::mutate(
      amd_clusts = as.integer(map_rank[as.character(amd_raw)]),
      amd_clusts = factor(amd_clusts, levels = seq_len(k), labels = paste0("CTS", seq_len(k)))
    )
}

revision_cts_occupancy <- function(cts_df) {
  cts_df %>%
    dplyr::filter(is.finite(age_ce), age_ce > 0) %>%
    dplyr::count(lake, amd_clusts, name = "n_samples") %>%
    dplyr::group_by(lake) %>%
    dplyr::mutate(prop = n_samples / sum(n_samples)) %>%
    dplyr::ungroup()
}

revision_write_csv <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(x, path)
  message("Wrote ", path)
  invisible(path)
}
