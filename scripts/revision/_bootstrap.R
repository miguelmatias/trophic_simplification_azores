# Shared bootstrap for revision A-scripts
revision_bootstrap <- function() {
  suppressPackageStartupMessages({
    library(tidyverse)
    library(vegan)
    library(mgcv)
    library(e1071)
    library(gratia)
  })

  root <- if (file.exists("data/clean_source_data_files.RData")) {
    "."
  } else if (file.exists("../../data/clean_source_data_files.RData")) {
    "../.."
  } else {
    stop("Run from repository root")
  }
  setwd(root)

  source("functions/source_custom_functions.R")
  source("functions/revision_helpers.R")
  revision_ensure_dirs()

  shared_path <- "outputs/revision/revision_shared.rds"
  if (!file.exists(shared_path)) {
    message("Shared objects missing; running setup...")
    source("scripts/revision/00_setup_revision_data.R")
  }
  readRDS(shared_path)
}
