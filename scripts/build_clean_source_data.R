#!/usr/bin/env Rscript
# Merge intermediate workspaces into data/clean_source_data_files.RData.
#
# Retains (from an optional full merged workspace) the same objects as the old
# make_for_publication.R subset — needed to recompute ordination in main_script.Rmd:
#   ls_df_diat_wide_codes, ls_df_codes_diat, ls_df_chiro_wide_codes,
#   df_chiro_codes, df_fgroups
#
# Also merges in gzip-compressed RData if present at repo root
# `data_files_comm_amd_globi.R` or `data/data_files_comm_amd_globi.R` (community AMD /
# GloBI score lists), then any objects already in data/clean_source_data_files.RData.
#
# Usage (from repository root):
#   Rscript scripts/build_clean_source_data.R
#   Rscript scripts/build_clean_source_data.R --merged path/to/data_files_guilds_globi_merged.R

args <- commandArgs(trailingOnly = TRUE)
merged_path <- NULL
publication <- FALSE
i <- 1L
while (i <= length(args)) {
  if (identical(args[[i]], "--merged") && i < length(args)) {
    merged_path <- args[[i + 1L]]
    i <- i + 2L
  } else if (identical(args[[i]], "--publication")) {
    publication <- TRUE
    i <- i + 1L
  } else {
    i <- i + 1L
  }
}

keep_from_merged <- c(
  "ls_df_diat_wide_codes",
  "ls_df_codes_diat",
  "ls_df_chiro_wide_codes",
  "df_chiro_codes",
  "df_fgroups"
)

prune_publication_inputs <- function(dest_env) {
  # Keep only columns needed to recompute ordinations in main_script.Rmd.
  #
  # Strategy:
  # - Wide species matrices: keep just identifiers needed for grouping/joining + taxa columns.
  # - Codes tables: keep only identifiers + age + (for diatoms) dec_depth used for comma fix.
  # - Functional groups: keep IDs + age_ce + group columns (drop age_bp).

  prune_wide <- function(df, id_keep) {
    id_keep <- intersect(id_keep, names(df))
    taxa_cols <- setdiff(names(df), names(df)[names(df) %in% c("lake", "core_id", "sec_cor", "samp_id", "samp_dep", "dec_depth", "core_depth_id", "depth_top", "depth_bot", "age_bp", "age_ce")])
    keep <- unique(c(id_keep, taxa_cols))
    df[, keep, drop = FALSE]
  }

  prune_codes <- function(df, keep) {
    keep <- intersect(keep, names(df))
    df[, keep, drop = FALSE]
  }

  if (exists("ls_df_diat_wide_codes", envir = dest_env, inherits = FALSE)) {
    x <- get("ls_df_diat_wide_codes", envir = dest_env)
    x <- lapply(x, prune_wide, id_keep = c("lake", "core_depth_id"))
    assign("ls_df_diat_wide_codes", x, envir = dest_env)
  }

  if (exists("ls_df_chiro_wide_codes", envir = dest_env, inherits = FALSE)) {
    x <- get("ls_df_chiro_wide_codes", envir = dest_env)
    x <- lapply(x, prune_wide, id_keep = c("lake", "core_depth_id", "depth_top"))
    assign("ls_df_chiro_wide_codes", x, envir = dest_env)
  }

  if (exists("ls_df_codes_diat", envir = dest_env, inherits = FALSE)) {
    x <- get("ls_df_codes_diat", envir = dest_env)
    x <- lapply(x, prune_codes, keep = c("lake", "core_depth_id", "dec_depth", "age_ce"))
    assign("ls_df_codes_diat", x, envir = dest_env)
  }

  if (exists("df_chiro_codes", envir = dest_env, inherits = FALSE)) {
    x <- get("df_chiro_codes", envir = dest_env)
    x <- lapply(x, prune_codes, keep = c("lake", "core_depth_id", "depth_top", "dec_depth", "age_ce"))
    assign("df_chiro_codes", x, envir = dest_env)
  }

  if (exists("df_fgroups", envir = dest_env, inherits = FALSE)) {
    df <- get("df_fgroups", envir = dest_env)
    group_cols <- c("high_profile", "low_profile", "motile", "euplanktonic", "algivore", "detritivore", "plantivore", "predator")
    keep <- intersect(c("lake", "core_depth_id", "age_ce", group_cols), names(df))
    assign("df_fgroups", df[, keep, drop = FALSE], envir = dest_env)
  }
}

cmdargs_trailingFalse <- commandArgs(trailingOnly = FALSE)
m <- grep("^--file=", cmdargs_trailingFalse, value = TRUE)
if (length(m)) {
  script_path <- normalizePath(sub("^--file=", "", m[[1L]]), mustWork = TRUE)
  root <- dirname(dirname(script_path))
} else {
  root <- getwd()
}
setwd(root)
message("Repository root: ", root)

outfile <- file.path(root, "data", "clean_source_data_files.RData")
in_clean <- file.path(root, "data", "clean_source_data_files.RData")

dest <- new.env(parent = emptyenv())

aux_names <- c("df_fgroups_glob_div", "cons_glob_div_plt_data", "prod_glob_div_plt_data")

if (publication) {
  if (is.null(merged_path)) {
    stop("For --publication builds you must provide --merged <path>.", call. = FALSE)
  }
  merged_path <- normalizePath(merged_path, mustWork = TRUE)
  e_m <- new.env(parent = emptyenv())
  load(merged_path, envir = e_m)
  missing <- setdiff(keep_from_merged, ls(e_m, all.names = FALSE))
  if (length(missing) > 0) {
    stop(
      "Merged file is missing required objects: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  for (nm in keep_from_merged) {
    assign(nm, get(nm, envir = e_m, inherits = FALSE), envir = dest)
  }
  prune_publication_inputs(dest)
  message("Copied publication subset from: ", merged_path)

  if (file.exists(in_clean)) {
    e_aux <- new.env(parent = emptyenv())
    load(in_clean, envir = e_aux)
    for (nm in aux_names) {
      if (exists(nm, envir = e_aux, inherits = FALSE)) {
        assign(nm, get(nm, envir = e_aux, inherits = FALSE), envir = dest)
      }
    }
    message("Copied auxiliary objects (if present) from: ", in_clean)
  }
} else {
  if (file.exists(in_clean)) {
    load(in_clean, envir = dest)
  }

  comm_candidates <- file.path(
    root,
    c(
      "data_files_comm_amd_globi.R",
      "data/data_files_comm_amd_globi.R"
    )
  )
  comm_path <- comm_candidates[file.exists(comm_candidates)][1]
  if (!is.na(comm_path)) {
    e_comm <- new.env(parent = emptyenv())
    con <- gzfile(comm_path, "rb")
    load(con, envir = e_comm)
    close(con)
    for (nm in ls(e_comm, all.names = FALSE)) {
      assign(nm, get(nm, envir = e_comm, inherits = FALSE), envir = dest)
    }
    message("Merged objects from gzip workspace: ", comm_path)
  }

  if (!is.null(merged_path)) {
    merged_path <- normalizePath(merged_path, mustWork = TRUE)
    e_m <- new.env(parent = emptyenv())
    load(merged_path, envir = e_m)
    missing <- setdiff(keep_from_merged, ls(e_m, all.names = FALSE))
    if (length(missing) > 0) {
      stop(
        "Merged file is missing required objects: ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
    }
    for (nm in keep_from_merged) {
      assign(nm, get(nm, envir = e_m, inherits = FALSE), envir = dest)
    }
    message("Copied publication subset from: ", merged_path)
  }
}

nms <- ls(dest, all.names = FALSE)
ordered <- unique(c(
  intersect(keep_from_merged, nms),
  setdiff(nms, keep_from_merged)
))

dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
save(list = ordered, envir = dest, file = outfile)
message(
  "Wrote ", outfile, " (", length(ordered), " objects): ",
  paste(ordered, collapse = ", ")
)
