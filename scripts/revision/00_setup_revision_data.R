# Setup shared revision objects for A1–A9
# Run from repo root: Rscript scripts/revision/00_setup_revision_data.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(mgcv)
  library(e1071)
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

message("Loading raw data...")
load("data/clean_source_data_files.RData")

norm_abund_guilds <- read.csv("data/table_abundance_guilds.csv", stringsAsFactors = FALSE)
nao_hernandez <- read.csv("data/table_nao_hernandez.csv", check.names = FALSE, stringsAsFactors = FALSE)
pollen <- read.csv("data/table_pollen_azores.csv", stringsAsFactors = FALSE)
lake_meta <- revision_load_lake_meta()

message("Building vegetation series (HGAM)...")
veg <- revision_build_vegetation_series(pollen)

message("Building 30-yr joined table...")
joined <- revision_build_joined_30yr(norm_abund_guilds, nao_hernandez, veg$veg_yearly)

message("Building arboreal + indicator onsets...")
arboreal <- revision_arboreal_series(pollen)
indicators <- revision_indicator_onsets(pollen)

lake_island <- lake_meta %>%
  dplyr::select(lake, island, alt, area, zmax)

shared <- list(
  norm_abund_guilds = norm_abund_guilds,
  lake_meta = lake_meta,
  lake_island = lake_island,
  veg_yearly = veg$veg_yearly,
  veg_dca = veg$dca,
  bio_lake_30 = joined$bio_lake_30,
  bio_reg_30 = joined$bio_reg_30,
  nao_30 = joined$nao_30,
  veg_30 = joined$veg_30,
  joined_df_30yr = joined$joined_df_30yr,
  arboreal = arboreal,
  indicators = indicators,
  ls_df_diat_wide_codes = ls_df_diat_wide_codes,
  ls_df_codes_diat = ls_df_codes_diat,
  ls_df_chiro_wide_codes = ls_df_chiro_wide_codes,
  df_chiro_codes = df_chiro_codes,
  df_fgroups = df_fgroups,
  cons_glob_div_plt_data = if (exists("cons_glob_div_plt_data")) cons_glob_div_plt_data else NULL,
  prod_glob_div_plt_data = if (exists("prod_glob_div_plt_data")) prod_glob_div_plt_data else NULL
)

outdir <- "outputs/revision"
saveRDS(shared, file.path(outdir, "revision_shared.rds"))
revision_write_csv(joined$joined_df_30yr, file.path(outdir, "joined_df_30yr.csv"))
revision_write_csv(veg$veg_yearly, file.path(outdir, "vegetation_yearly_smooth.csv"))
revision_write_csv(arboreal, file.path(outdir, "arboreal_pct_by_site.csv"))
revision_write_csv(indicators, file.path(outdir, "indicator_onsets.csv"))

message("Setup complete: outputs/revision/revision_shared.rds")
