# Consolidate uploaded tephra events + references into
# data/revision/tephra_events.csv (canonical A6 input).
suppressPackageStartupMessages(library(dplyr))

read_utf8_csv <- function(path) {
  if (requireNamespace("readr", quietly = TRUE)) {
    df <- readr::read_csv(
      path,
      show_col_types = FALSE,
      locale = readr::locale(encoding = "UTF-8")
    )
    return(as.data.frame(df, stringsAsFactors = FALSE))
  }
  raw <- readBin(path, what = "raw", n = file.info(path)$size)
  if (length(raw) >= 3 && raw[1] == as.raw(0xef) && raw[2] == as.raw(0xbb) && raw[3] == as.raw(0xbf)) {
    raw <- raw[-(1:3)]
  }
  txt <- rawToChar(raw)
  Encoding(txt) <- "UTF-8"
  utils::read.csv(text = txt, stringsAsFactors = FALSE, check.names = FALSE)
}

events <- read_utf8_csv("data/revision/tephra_events_raw.csv")
refs <- read_utf8_csv("data/revision/tephra_references.csv")

# Drop regional-only provenance rows (lake like "none (regional)"); keep lake-specific inventory
events <- events %>%
  filter(!grepl("^none", lake, ignore.case = TRUE))

primary_ref <- trimws(sub(";.*$", "", events$reference))

refs_use <- refs %>%
  filter(!grepl("^Secondary", short_ref)) %>%
  transmute(
    primary_ref = short_ref,
    doi = doi,
    full_citation = full_citation,
    ref_checked = checked
  )

out <- events %>%
  mutate(
    lake = as.character(lake),
    lake = iconv(lake, from = "", to = "UTF-8", sub = ""),
    lake = trimws(lake),
    lake = gsub("Caldeir.?o", "Caldeirao", lake, perl = TRUE),
    lake = gsub("Emp\\.\\s*Norte", "Empadadas Norte", lake),
    lake = gsub("^Emp Norte$", "Empadadas Norte", lake),
    age_ce = as.numeric(age_ce),
    primary_ref = primary_ref,
    # Lake-specific SUPPORTED layers only (drop proxy-inferred TENTATIVE / NOT SUPPORTED)
    include_sensitivity =
      grepl("^SUPPORTED", status, ignore.case = TRUE) &
      is.finite(age_ce)
  ) %>%
  left_join(refs_use, by = "primary_ref") %>%
  select(
    lake, island, age_ce, age_as_published, event, status,
    evidence_in_lake, source, reference, primary_ref, doi, full_citation,
    ref_checked, notes, include_sensitivity
  ) %>%
  arrange(lake, age_ce)

# write with write.csv; ensure UTF-8
con <- file("data/revision/tephra_events.csv", open = "w", encoding = "UTF-8")
on.exit(close(con), add = TRUE)
write.csv(out, con, row.names = FALSE, na = "")

message("n rows: ", nrow(out))
message(
  "n lakes (excl regional): ",
  length(unique(out$lake[!grepl("^none", out$lake)]))
)
message("lakes: ", paste(sort(unique(out$lake)), collapse = ", "))
message("include_sensitivity TRUE: ", sum(out$include_sensitivity, na.rm = TRUE))
message("Unmatched primary_ref (no DOI):")
print(sort(unique(out$primary_ref[is.na(out$doi) | out$doi == ""])))
message("Sensitivity rows:")
sens <- out %>% filter(include_sensitivity) %>% select(lake, age_ce, status, event)
print(as.data.frame(sens), row.names = FALSE)
