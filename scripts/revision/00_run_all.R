# Run all revision analyses A1–A9
# Usage (from repo root): Rscript scripts/revision/00_run_all.R

scripts <- c(
  "scripts/revision/00_setup_revision_data.R",
  "scripts/revision/A04_nao_vegetation_correlation.R",
  "scripts/revision/A07_split_cts5_cts6.R",
  "scripts/revision/A09_onset_vs_stocking.R",
  "scripts/revision/A01_rarefy_richness.R",
  "scripts/revision/A02_resolution_matched.R",
  "scripts/revision/A03_common_period_coverage.R",
  "scripts/revision/A05_human_vs_climate_vegetation.R",
  "scripts/revision/A06_volcanism_sensitivity.R",
  "scripts/revision/A08_morphometry_cts.R",
  "scripts/revision/A10_fig3e_rarefied_cts.R",
  "scripts/revision/A11_alt_figure4_tephra.R"
)

# Avoid re-running setup twice when individual scripts bootstrap
run_one <- function(path) {
  message("\n############################\n# Running ", path, "\n############################")
  # Each A-script sources bootstrap; setup script is standalone
  if (basename(path) == "00_setup_revision_data.R") {
    source(path, local = new.env())
  } else {
    # Run in fresh process-like local env but shared search path
    env <- new.env(parent = globalenv())
    sys.source(path, envir = env)
  }
  message("Finished ", path)
}

# Prefer sequential sourcing in current session for shared package state
if (!file.exists("data/clean_source_data_files.RData")) {
  stop("Run from repository root")
}

# Setup first
source("scripts/revision/00_setup_revision_data.R")

# Then analyses (skip setup entry)
for (s in scripts[-1]) {
  message("\n############################\n# Running ", s, "\n############################")
  source(s)
}

message("\nAll revision analyses finished. See outputs/revision/")
