## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Set-up chunk
knitr::opts_chunk$set(
  echo = FALSE,
  warning = FALSE,
  fig.pos = "H",
  fig.width = 12
)


## ----sources, include=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------

# -------------------------------
# 🧰 General-Purpose & Data Handling
# -------------------------------
library(tidyverse)      # Loads: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats
library(zoo)            # Time-series objects and interpolation
library(parallel)       # Base R parallel processing support
library(furrr)          # Parallel purrr (requires future backend)
library(sf)             # Spatial vector data handling

# -------------------------------
# 🖍 Aesthetic & Utility Enhancements
# -------------------------------
library(scales)         # Useful for custom ggplot2 axes
library(png)            # For reading PNG images (only if you're embedding raster images)
library(ggtext)         # Rich text/markdown/emoji in ggplot2 elements
library(Hmisc)          # For `capitalize()` or advanced summary stats
library(scico)          # Color palettes (colorblind-safe, perceptually uniform)
library(conflicted)     # Safer namespace conflict management

# -------------------------------
# 📊 Visualization & Plot Layout
# -------------------------------
# ggplot2 already loaded via tidyverse, but other extensions aren't:
library(ggrepel)        # Smart label positioning
library(ggpmisc)        # Annotations and statistical labels
library(viridis)        # Better color palettes
library(cowplot)        # Clean multi-panel figures
library(patchwork)      # Combine ggplots with intuitive syntax
library(gridExtra)      # Arrange base or ggplots in grids
library(grid)           # Base graphics system (already available, no need to load)
library(ggnewscale)     # Multiple color/fill scales
library(ggridges)       # Ridgeline plots
library(ggh4x)          # Extended faceting and axes
library(ggspatial)      # North arrows, scale bars, etc.
library(ggordiplots)    # Ordination diagrams (RDA, CCA, NMDS)

# -------------------------------
# 📈 Statistical Analysis & Modeling
# -------------------------------
library(vegan)          # Ecological ordination (NMDS, CCA, RDA), diversity metrics
library(codyn)          # Temporal beta-diversity and dynamics
library(broom)          # Tidy model outputs
library(mgcv)           # Generalized Additive Models
library(scam)           # Shape-constrained additive models
library(gratia)         # Diagnostics and visualization for `mgcv::gam()`

# -------------------------------
# 🤖 Machine Learning & Interpretation
# -------------------------------
library(caret)          # Framework for ML model training/testing
library(e1071)          # SVMs, Naive Bayes, etc. (caret dependency)
library(evtree)         # Evolutionary trees
library(rpart)          # Classification/regression trees
library(randomForest)   # Random Forests
library(rpart.plot)     # Visualize rpart trees
library(earth)          # MARS regression
library(DALEX)          # Model explanation tools
library(DALEXtra)       # Extra methods for xgboost, caret, etc.
library(irr)            # Inter-rater agreement (Cohen's Kappa, etc.)
library(ape)            # Phylogenetic tools (tree manipulation, distances)

# -------------------------------
# 📂 Custom Functions
# -------------------------------
source("functions/source_custom_functions.R")  # Custom paleoecological tools

# -------------------------------
# 🚨 Solve Function Conflicts
# -------------------------------
conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("summarise", "dplyr")
conflict_prefer("summarize", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("explain", "DALEX")
conflict_prefer("annotate", "ggplot2")


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Load local pre-saved RData files (contains trophic data and clustering outputs)
# ------------------------------------------------------------------------------
load("data/data_files.R")
load("data/data_files_comm_amd_globi.R")
load("data/data_files_guilds_globi_merged.R")

# Define functional groups (used throughout trophic structure analyses)
alternative_functional_groups_guilds_globi_merged <- c(
  "high_profile", "low_profile", "motile", "euplanktonic",
  "algivore", "detritivore", "plantivore", "predator"
)

# ------------------------------------------------------------------------------
# Load lake metadata
# ------------------------------------------------------------------------------
lake_metadata <- read.csv(file = "data/table_lake_metadata.csv")  # lake coordinates, area, etc.

# ------------------------------------------------------------------------------
# Load guild relative abundances and diversity tables
# ------------------------------------------------------------------------------
norm_abund_guilds <- read.table("data/table_abundance_guilds.csv",header = TRUE, sep = ",", dec = ".")
div_guilds <- read.table("data/table_diversity_guilds.csv", header = TRUE, sep = ",", dec = ".")

# ------------------------------------------------------------------------------
# Load climate and environmental reconstructions
# ------------------------------------------------------------------------------
nao_hernandez <- read.table("data/table_nao_hernandez.csv", header = TRUE, sep = ",", dec = ".")
nhst_buntgen <- read.table("data/table_nhst_buntgen.csv", header = TRUE, sep = ",", dec = ".")

# ------------------------------------------------------------------------------
# Load pollen group data
# ------------------------------------------------------------------------------
pollen <- read.table("data/table_pollen_azores.csv", header = TRUE, sep = ",", dec = ".")

## ----fig.asp = 1, fig.width = 12----------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Extract PCoA vectors (1st and 2nd axes) for each lake
# ------------------------------------------------------------------------------

ls_df_pcoa_scores_diatoms_hellinger <- ls_df_diat_wide_codes %>%
  # Remove metadata columns
  map(., ~ .x %>% select(-c(lake:core_depth_id))) %>%
  # Replace NAs with zeros
  map(., ~ .x %>% replace(is.na(.), 0)) %>%
  # Remove species with total abundance = 0 across all samples
  map(., ~ .x %>% select_if(list(~ sum(.x) > 0))) %>%
  # Apply Hellinger transformation
  map(., ~ .x %>% decostand(method = "hellinger")) %>%
  # Compute Bray–Curtis dissimilarity
  map(., ~ .x %>% vegdist(method = "bray", na.rm = TRUE)) %>%
  # Run PCoA
  map(., ~ .x %>% pcoa()) %>%
  # Extract eigenvectors
  map(., function(x) x$vectors) %>%
  # Convert first two axes into a data frame
  map(., function(x) as.data.frame(unlist(x[, 1:2])))

# ------------------------------------------------------------------------------
# Bind metadata and PCoA scores for each lake
# ------------------------------------------------------------------------------

ls_df_all_scores_diatoms_amd_globi <- NULL

for (i in 1:length(ls_df_pcoa_scores_diatoms_hellinger)) {
  ls_df_all_scores_diatoms_amd_globi[[i]] <- bind_cols(
    ls_df_codes_diat[[i]],
    ls_df_pcoa_scores_diatoms_hellinger[[i]]
  )
}

# ------------------------------------------------------------------------------
# Fix decimal comma issue in depth field
# ------------------------------------------------------------------------------

ls_df_all_scores_diatoms_amd_globi <- ls_df_all_scores_diatoms_amd_globi %>%
  map(., ~ .x %>% mutate(dec_depth = as.numeric(str_replace(dec_depth, ",", "."))))

# ------------------------------------------------------------------------------
# Calculate abundance, richness, and Simpson index for each sample
# ------------------------------------------------------------------------------

ls_df_diatoms_wide_abund <- ls_df_diat_wide_codes %>%
  # Select species-only columns
  map(., ~ .x %>% select(8:last_col())) %>%
  # Add total abundance column
  map(., ~ .x %>% mutate(lake_abundance = rowSums(.x))) %>%
  # Add richness (count of present species)
  map(., ~ .x %>% mutate(local_richness = apply(.x[8:(ncol(.x) - 1)] > 0, 1, sum))) %>%
  # Add Simpson diversity index
  map(., ~ .x %>% mutate(simpson = apply(.x[8:(ncol(.x) - 1)], 1, function(x) abdiv::simpson(x)))) %>%
  # Keep only calculated metrics
  map(., ~ .x %>% select(lake_abundance, local_richness, simpson))

# ------------------------------------------------------------------------------
# Bind PCoA scores with richness and abundance metrics
# ------------------------------------------------------------------------------

for (i in 1:length(ls_df_diatoms_wide_abund)) {
  ls_df_all_scores_diatoms_amd_globi[[i]] <- bind_cols(
    ls_df_all_scores_diatoms_amd_globi[[i]],
    ls_df_diatoms_wide_abund[[i]]
  )
}

# ------------------------------------------------------------------------------
# Extract DCA vectors (1st and 2nd axes) for each lake
# ------------------------------------------------------------------------------

ls_dca_diatoms <- ls_df_diat_wide_codes %>%
  # Remove metadata columns
  map(., ~ .x %>% select(-c(lake:core_depth_id))) %>%
  # Remove species with zero abundance across all samples
  map(., ~ .x %>% select_if(list(~ sum(.x) > 0))) %>%
  # Apply Hellinger transformation
  map(., ~ .x %>% decostand(method = "hellinger")) %>%
  # Perform DCA
  map(., ~ .x %>% decorana()) %>%
  # Extract DCA scores
  map(., function(x) scores(x)) %>%
  # Convert first two axes into a data frame
  map(., function(x) as.data.frame(unlist(x[, 1:2])))

# ------------------------------------------------------------------------------
# Bind DCA scores to existing PCoA + richness/abundance data
# ------------------------------------------------------------------------------

for (i in 1:length(ls_df_all_scores_diatoms_amd_globi)) {
  ls_df_all_scores_diatoms_amd_globi[[i]] <- bind_cols(
    ls_df_all_scores_diatoms_amd_globi[[i]],
    ls_dca_diatoms[[i]]
  )
}

# ------------------------------------------------------------------------------
# Create additional list with DCA scores + metadata
# ------------------------------------------------------------------------------

ls_dca_diatoms_2 <- ls_df_diat_wide_codes %>%
  map(., ~ .x %>% select(-c(lake:core_depth_id))) %>%
  map(., ~ .x %>% select_if(list(~ sum(.x) > 0))) %>%
  map(., ~ .x %>% decostand(method = "hellinger")) %>%
  map(., ~ .x %>% decorana()) %>%
  map(., function(x) scores(x)) %>%
  map(., function(x) as.data.frame(unlist(x[, 1:2])))

for (i in 1:length(ls_dca_diatoms_2)) {
  ls_dca_diatoms_2[[i]] <- bind_cols(
    ls_df_codes_diat[[i]],
    ls_dca_diatoms_2[[i]]
  )
}



## ----fig.asp = 1, fig.width = 12----------------------------------------------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Extract PCoA vectors (1st and 2nd axes) for each lake based on chironomid data
# ------------------------------------------------------------------------------

ls_df_pcoa_scores_chiro_hellinger <- ls_df_chiro_wide_codes %>%
  # Remove metadata columns
  map(., ~ .x %>% select(-c(lake:depth_top))) %>%
  # Remove species columns with no occurrences
  map(., ~ .x %>% select_if(list(~ sum(.x) > 0))) %>%
  # Apply Hellinger transformation
  map(., ~ .x %>% decostand(method = "hellinger")) %>%
  # Compute Bray–Curtis dissimilarity matrix
  map(., ~ .x %>% vegdist(method = "bray", na.rm = TRUE)) %>%
  # Perform PCoA
  map(., ~ .x %>% pcoa()) %>%
  # Extract eigenvectors
  map(., function(x) x$vectors) %>%
  # Keep first two axes and convert to data frame
  map(., function(x) as.data.frame(unlist(x[, 1:2])))

# ------------------------------------------------------------------------------
# Combine metadata with PCoA scores
# ------------------------------------------------------------------------------

ls_df_all_scores_chiro_amd_globi <- NULL

for (i in 1:length(ls_df_pcoa_scores_chiro_hellinger)) {
  ls_df_all_scores_chiro_amd_globi[[i]] <- bind_cols(
    df_chiro_codes[[i]],
    ls_df_pcoa_scores_chiro_hellinger[[i]]
  )
}

# ------------------------------------------------------------------------------
# Assign lake names to list entries
# ------------------------------------------------------------------------------

names(ls_df_all_scores_chiro_amd_globi) <- unlist(get_names_list(ls_df_all_scores_chiro_amd_globi))

# ------------------------------------------------------------------------------
# Calculate chironomid abundance, richness, and Simpson index
# ------------------------------------------------------------------------------

lake_abund_chiro <- ls_df_chiro_wide_codes %>%
  # Keep only species columns
  map(., ~ .x %>% select(8:last_col())) %>%
  # Compute total abundance per sample
  map(., ~ .x %>% mutate(chiro_lake_abundance = rowSums(.x))) %>%
  # Compute local richness (number of present species)
  map(., ~ .x %>% mutate(chiro_local_richness = apply(.x[8:(ncol(.x) - 1)] > 0, 1, sum))) %>%
  # Compute Simpson diversity index
  map(., ~ .x %>% mutate(chiro_simpson = apply(.x[8:(ncol(.x) - 1)], 1, function(x) abdiv::simpson(x)))) %>%
  # Retain calculated metrics only
  map(., ~ .x %>% select(chiro_lake_abundance, chiro_local_richness, chiro_simpson))

# ------------------------------------------------------------------------------
# Bind abundance and diversity metrics to main list
# ------------------------------------------------------------------------------

for (i in 1:length(ls_df_all_scores_chiro_amd_globi)) {
  ls_df_all_scores_chiro_amd_globi[[i]] <- bind_cols(
    ls_df_all_scores_chiro_amd_globi[[i]],
    lake_abund_chiro[[i]]
  )
}

# ------------------------------------------------------------------------------
# Extract DCA vectors (1st and 2nd axes) from chironomid data
# ------------------------------------------------------------------------------

ls_dca_chiro <- ls_df_chiro_wide_codes %>%
  # Remove metadata columns
  map(., ~ .x %>% select(-c(lake:depth_top))) %>%
  # Remove species columns with no data
  map(., ~ .x %>% select_if(list(~ sum(.x) > 0))) %>%
  # Apply Hellinger transformation
  map(., ~ .x %>% decostand(method = "hellinger")) %>%
  # Perform DCA
  map(., ~ .x %>% decorana()) %>%
  # Extract site scores
  map(., function(x) scores(x)) %>%
  # Keep first two axes and convert to data frame
  map(., function(x) as.data.frame(unlist(x[, 1:2])))

# ------------------------------------------------------------------------------
# Bind DCA scores to metadata + PCoA + diversity info
# ------------------------------------------------------------------------------

for (i in 1:length(ls_df_all_scores_chiro_amd_globi)) {
  ls_df_all_scores_chiro_amd_globi[[i]] <- bind_cols(
    ls_df_all_scores_chiro_amd_globi[[i]],
    ls_dca_chiro[[i]]
  )
}

# ------------------------------------------------------------------------------
# Create second DCA list with metadata attached
# ------------------------------------------------------------------------------

ls_dca_chiro_2 <- ls_df_chiro_wide_codes %>%
  map(., ~ .x %>% select(-c(lake:core_depth_id))) %>%
  map(., ~ .x %>% select_if(list(~ sum(.x) > 0))) %>%
  map(., ~ .x %>% decostand(method = "hellinger")) %>%
  map(., ~ .x %>% decorana()) %>%
  map(., function(x) scores(x)) %>%
  map(., function(x) as.data.frame(unlist(x[, 1:2])))

for (i in 1:length(ls_dca_chiro_2)) {
  ls_dca_chiro_2[[i]] <- bind_cols(
    df_chiro_codes[[i]],
    ls_dca_chiro_2[[i]]
  )
}




## ----fig.asp = 1, fig.width = 12----------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Clean initial data and remove NA rows
# ------------------------------------------------------------------------------

df_fgroups <- df_fgroups[!is.na(rowSums(df_fgroups[, -c(1:4)])), ]

# Ungroup any groupings that might exist in df_fgroups
df_fg <- df_fgroups %>% ungroup()

# Extract identifying metadata (lake, core_id, depth, age_ce)
df_fg_codes <- df_fg %>%
  ungroup() %>%
  select(lake:age_ce) %>%
  distinct()

# Split the metadata by lake into a list
ls_df_fgroups_wide_codes <- df_fg_codes %>%
  group_split(lake)

# Assign names to list based on lake
names(ls_df_fgroups_wide_codes) <- unlist(get_names_list(ls_df_fgroups_wide_codes))

# ------------------------------------------------------------------------------
# Run DCA on each lake's functional group composition
# ------------------------------------------------------------------------------

ls_dca_fgroups <- df_fg %>%
  group_split(lake) %>%
  # Remove metadata columns
  map(., ~ .x %>% dplyr::select(-c(lake:age_ce))) %>%
  # Remove empty rows (samples with no groups)
  map(., ~ .x %>% filter(rowSums(.) != 0, na.rm = TRUE)) %>%
  # Remove empty group columns (groups not present)
  map(., ~ .x %>% select_if(list(~ sum(.x) > 0))) %>%
  # Apply Hellinger transformation
  map(., ~ .x %>% decostand(method = "hellinger")) %>%
  # Perform DCA ordination
  map(., ~ .x %>% decorana()) %>%
  # Extract sample scores
  map(., function(x) scores(x)) %>%
  # Convert first two axes to data frame
  map(., function(x) as.data.frame(unlist(x[, 1:2])))

# Also create a list of raw data split by lake (used later)
ls_dca_fgroups_EWS <- df_fg %>% group_split(lake)

# ------------------------------------------------------------------------------
# Combine metadata and DCA results
# ------------------------------------------------------------------------------

# Combine all DCA outputs into one data frame and bind with metadata
df_dca_fgroups_codes_wide <- ls_dca_fgroups %>%
  bind_rows() %>%
  bind_cols(df_fg_codes, .)

# Split DCA results with metadata by lake
ls_df_dca_fgroups_codes_wide <- df_dca_fgroups_codes_wide %>%
  group_split(lake)

# Assign names to list
names(ls_df_dca_fgroups_codes_wide) <- unlist(get_names_list(ls_df_dca_fgroups_codes_wide))

# ------------------------------------------------------------------------------
# Prepare input for PCoA analyses
# ------------------------------------------------------------------------------

# Split raw abundance data by lake
ls_df_fg <- df_fg %>%
  group_split(lake)

# Assign names to list
names(ls_df_fg) <- unlist(get_names_list(ls_df_fg))

# ------------------------------------------------------------------------------
# Compute PCoA for each lake using Bray–Curtis distances
# ------------------------------------------------------------------------------

ls_df_pcoa_scores_fgroups <- ls_df_fg %>%
  # Remove metadata columns
  map(., ~ .x %>% select(-c(lake:age_ce))) %>%
  # Replace NAs with 0s
  map(., ~ .x %>% replace(is.na(.), 0)) %>%
  # Remove empty group columns
  map(., ~ .x %>% select_if(list(~ sum(.x) > 0))) %>%
  # Optional: Apply Hellinger transformation (common before Bray-Curtis, but not mandatory here)
  map(., ~ .x %>% decostand(method = "hellinger")) %>%
  # Calculate Bray–Curtis distance
  map(., ~ .x %>% vegdist(method = "bray", na.rm = TRUE)) %>%
  # Perform PCoA
  map(., ~ .x %>% pcoa()) %>%
  # Extract sample scores
  map(., function(x) x$vectors) %>%
  # Keep first two PCoA axes and convert to data frame
  map(., function(x) as.data.frame(unlist(x[, 1:2])))

# ------------------------------------------------------------------------------
# Bind metadata and PCoA scores
# ------------------------------------------------------------------------------

for (i in 1:length(ls_df_pcoa_scores_fgroups)) {
  ls_df_pcoa_scores_fgroups[[i]] <- bind_cols(
    ls_df_dca_fgroups_codes_wide[[i]],
    ls_df_pcoa_scores_fgroups[[i]]
  )
}



## ----fig.asp = 1, fig.width = 12----------------------------------------------------------------------------------------------------------------------------------------------------------------

# --------------------------------------------
# Apply axis correction to diatom DCA scores
# --------------------------------------------
b_ls_dca_diatoms <- 
  ls_dca_diatoms_2 %>%
  map(~ .x %>% arrange(age_ce)) %>%         # Sort by age
  map(~ .x %>% convert_DCA()) %>%           # Flip DCA if needed (Custom function in functions file)
  bind_rows() %>%
  select(lake, age_ce, DCA1, DCA2)

# Fix Prata manually (function failed or returned incorrect orientation)
b_ls_dca_diatoms[b_ls_dca_diatoms$lake == "Prata", "DCA1"] <- 
  b_ls_dca_diatoms[b_ls_dca_diatoms$lake == "Prata", "DCA1"] * -1

# --------------------------------------------
# Apply axis correction to chironomid DCA scores
# --------------------------------------------
b_ls_dca_chiro <- 
  ls_dca_chiro_2 %>%
  map(~ .x %>% arrange(age_ce)) %>%
  map(~ .x %>% convert_DCA()) %>%
  bind_rows() %>%
  select(lake, age_ce, DCA1, DCA2)

# --------------------------------------------
# Apply axis correction to functional group DCA scores
# --------------------------------------------
b_ls_dca_fgroup <- 
  ls_df_dca_fgroups_codes_wide %>%
  map(~ .x %>% arrange(age_ce)) %>%
  map(~ .x %>% convert_DCA()) %>%
  bind_rows() %>%
  select(lake, age_ce, DCA1, DCA2)

# --------------------------------------------
# Combine all DCA scores into a single long-format dataframe
# --------------------------------------------
df_dca_multi <- bind_rows(
  b_ls_dca_diatoms,
  b_ls_dca_chiro,
  b_ls_dca_fgroup,
  .id = "troph"  # 1 = diatoms, 2 = chiro, 3 = fgroups
) %>%
  mutate(group = as.factor(troph))  # Ensure group is treated as factor for plotting

# --------------------------------------------
# Split data by lake into a list for faceted plotting
# --------------------------------------------
ls_dca_multi <- df_dca_multi %>%
  group_split(lake)

# Assign lake names as list names
names(ls_dca_multi) <- unlist(get_names_list(ls_dca_multi))

# --------------------------------------------
# Create list of plots comparing diatoms vs chironomids by lake
# --------------------------------------------
ls_plot_df_dca_multi_all <- 
  ls_dca_multi %>%
  map(~ .x %>% filter(group %in% 1:2)) %>%  # Filter to diatoms (1) and chironomids (2)
  map2(
    ., names(ls_dca_multi),
    ~ ggplot(.x, aes(x = age_ce, y = DCA1)) +
      geom_smooth(
        aes(colour = group), 
        linewidth = 2, alpha = 0.5, method = "gam", show.legend = TRUE
      ) +
      # Optional: Add raw points if desired
      # geom_point(aes(color = group), size = 3) +
      scale_colour_viridis_d(
        name = "Groups",
        labels = c("Producers", "Consumers")
      ) +
      theme(plot.background = NULL) +
      labs(title = .y) +  # Lake name as plot title
      theme_bw() +
      theme(
        plot.title = element_text(face = "bold", size = 12),
        legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
        legend.justification = c(0, 1),
        legend.position.inside = c(0, 1),
        axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.minor = element_blank()
      )
  )




## ----fig.asp = 0.5, fig.width = 12--------------------------------------------------------------------------------------------------------------------------------------------------------------

# --------------------------------------------
# Combine list of functional group data frames into one table
# --------------------------------------------
tab_sums <- 
  ls_df_fg %>% 
  map(~.x %>% as.data.frame()) %>%  # Convert each list element to data frame
  bind_rows()  # Combine into one unified data frame

# --------------------------------------------
# Normalize abundance data by row (i.e., relative abundance)
# --------------------------------------------
norm_abund <- 
  bind_cols(
    select(tab_sums, lake:age_ce),  # Retain identifying metadata
    select(tab_sums, high_profile:euplanktonic) / rowSums(select(tab_sums, high_profile:euplanktonic)),  # Normalize functional groups
    select(tab_sums, algivore:last_col()) / rowSums(select(tab_sums, algivore:last_col()))  # Normalize trait groups
  ) %>% 
  replace(is.na(.), 0)  # Replace any resulting NAs with 0s

# --------------------------------------------
# Detect available CPU cores for parallel processing
# --------------------------------------------
numCores <- parallel::detectCores()

# --------------------------------------------
# 4. Run clustering optimization to identify best AMD cluster number
# --------------------------------------------
ls_opt_clusts <- 
  norm_abund %>%
  select(-c(lake:age_ce)) %>%  # Remove metadata columns
  PlotoptAMDclusters(., .iterations = 1000, .min_clusts = 3, .num_groups = 9)  # Evaluate from 3 to 9 clusters

# --------------------------------------------
# Identify optimal number of clusters via asymptote detection
# --------------------------------------------
x <- ls_opt_clusts$results$.clusters
y <- ls_opt_clusts$results$.maxpmean
dy <- diff(y)  # Rate of change in model performance
threshold <- 0.01  # Define a threshold for plateau detection
asymptote_point <- which(abs(dy) < threshold)[1]  # First point where improvement stabilizes

# Visual check for asymptote
plot(x, y, type = "l", main = "Identifying Asymptote")
points(x[asymptote_point], y[asymptote_point], col = "red", pch = 19, cex = 1.5)
abline(v = x[asymptote_point], col = "blue", lty = 2)
text(x[asymptote_point], y[asymptote_point], labels = paste0("x = ", x[asymptote_point]), pos = 4)

# Extract cluster number at asymptote
ls_opt_clusts <- as.vector(ls_opt_clusts$results[asymptote_point, 1])$.clusters

# --------------------------------------------
# Assign cluster identity to each row (sample)
# --------------------------------------------
opt_AMD_clusts_fgroups_paleo <- 
  norm_abund %>%
  select(-c(lake:age_ce)) %>%  # Remove metadata
  select_if(~ sum(.) > 0) %>%  # Drop zero-only columns
  getAMDclusters(.data = ., .opt_num_clusts = ls_opt_clusts, .iterations = 5000) %>% 
  dplyr::rename(amd_clusts = .clst)  # Rename cluster column

# Ensure cluster is a proper factor
opt_AMD_clusts_fgroups_paleo <- opt_AMD_clusts_fgroups_paleo %>%
  mutate(amd_clusts = factor(as.vector(unlist(.[, "amd_clusts"]))))

# --------------------------------------------
# Prepare ordination scores and metadata
# --------------------------------------------
glob_df_pcoa_scores_fgroups <- bind_rows(ls_df_pcoa_scores_fgroups)

# Subset ordination scores to match samples in norm_abund
sub_glob_df_pcoa_scores_fgroups <- 
  glob_df_pcoa_scores_fgroups %>% 
  filter(core_depth_id %in% norm_abund$core_depth_id)

# --------------------------------------------
# Merge metadata, ordination, and cluster assignment into one table
# --------------------------------------------
global_matrix_merged_paleo <- 
  bind_cols(
    select(norm_abund, lake:age_ce),
    select(sub_glob_df_pcoa_scores_fgroups, DCA1:DCA2),
    amd_clusts = factor(opt_AMD_clusts_fgroups_paleo[, 1])  # Add cluster assignment
  ) %>% 
  left_join(., 
            select(ungroup(df_fgroups_glob_div), core_depth_id:total_nspp_by_lake_core) %>% 
              distinct(),
            by = "core_depth_id")

# --------------------------------------------
# Create color mapping key based on cluster means
# --------------------------------------------
key <- 
  bind_cols(
    global_matrix_merged_paleo,
    select(norm_abund, c(high_profile:last_col()))
  ) %>% 
  group_by(amd_clusts) %>% 
  summarise(
    mean_euplanktonic = mean(euplanktonic, na.rm = TRUE),
    mean_total_nspp_by_lake_core = mean(total_nspp_by_lake_core, na.rm = TRUE)
  ) %>% 
  arrange(mean_total_nspp_by_lake_core) %>% 
  mutate(rank = row_number())  # Assign ranking based on richness

# --------------------------------------------
# Apply ranking to cluster IDs (relabel clusters by complexity)
# --------------------------------------------
cluster_mapping <- key %>% 
  select(amd_clusts, rank) %>% 
  deframe()  # Convert to named vector

# Replace cluster IDs with rank values
global_matrix_merged_paleo <- global_matrix_merged_paleo %>%
  mutate(amd_clusts = cluster_mapping[as.character(amd_clusts)]) %>%
  mutate(amd_clusts = as.factor(amd_clusts))

# --------------------------------------------
# [Decision made here] Merge two cluster levels
# --------------------------------------------
global_matrix_merged_paleo$amd_clusts[global_matrix_merged_paleo$amd_clusts == 6] <- 5

# --------------------------------------------
# Split final data frame by lake into a named list
# --------------------------------------------
ls_all_hclusts_complete_paleo <- 
  global_matrix_merged_paleo %>%
  group_split(lake)

# Assign names to list elements using helper function
names(ls_all_hclusts_complete_paleo) <- unlist(get_names_list(ls_all_hclusts_complete_paleo))




## ----fig.asp=1, fig.width=12--------------------------------------------------------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------
# Merge all lake-specific dataframes into one long-format table
# ----------------------------------------------
merged_data <- ls_all_hclusts_complete_paleo %>%
  map(~.x %>% as.data.frame()) %>%
  bind_rows() %>%
  filter(age_ce > 0) %>%  # Remove non-positive ages (pre-zero)
  select(lake, core_depth_id, amd_clusts, total_nspp_by_lake_core) %>%
  unique()

# ----------------------------------------------
# Perform ANOVA and Tukey post-hoc test
# ----------------------------------------------
anova_res <- aov(total_nspp_by_lake_core ~ amd_clusts, data = merged_data)
tukey_res <- TukeyHSD(anova_res)

# Generate compact letter display to indicate group differences
cld <- multcompView::multcompLetters4(anova_res, tukey_res)

# Extract significant letters and prepare for plotting
signif_letters <- data.frame(
  amd_clusts = names(cld$amd_clusts$Letters),
  letters = cld$amd_clusts$Letters
)

# Add y-position for each letter above boxplots
signif_letters <- signif_letters %>%
  left_join(
    merged_data %>%
      group_by(amd_clusts) %>%
      summarise(max_y = max(total_nspp_by_lake_core, na.rm = TRUE)),
    by = "amd_clusts"
  ) %>%
  mutate(y_position = max_y * 1.1)

# ----------------------------------------------
# Create boxplot comparing diversity across clusters (all lakes combined)
# ----------------------------------------------
box_plot_clusts <- ggplot(merged_data, aes(x = amd_clusts, y = total_nspp_by_lake_core)) +
  geom_boxplot(aes(fill = amd_clusts)) +
  geom_text(
    aes(x = amd_clusts, y = y_position, label = letters),
    inherit.aes = FALSE,
    data = signif_letters, size = 6
  ) +
  ggpubr::stat_compare_means(method = "anova", label.y = 9, cex = 4) +  # ANOVA p-value label
  scale_fill_viridis_d(name = "Community Trophic \nstructures", option = "viridis", direction = -1) +
  theme_minimal() +
  labs(x = "Community Trophic Structures", y = "Number of species") +
  theme(
    plot.title = element_text(size = 12),
    legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "bottom",
    axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Assuming your x-axis variable is called 'x_var' in fig_2_time_clusts data:
box_plot_clusts$data$amd_clusts <- factor(box_plot_clusts$data$amd_clusts, 
                                       levels = c(6,5,4,3,2,1))
box_plot_clusts_mod <- 
  box_plot_clusts + 
  scale_fill_viridis_d(direction = 1) +
  theme(legend.position = "none")

# Rebuild the plot to reflect the new x-axis order
box_plot_clusts_mod <- box_plot_clusts_mod +
  scale_x_discrete()  # Ensures x-axis reflects the new order

# Create the inset plot
box_plot_clusts_mod <- box_plot_clusts_mod +
  #theme(plot.background = element_rect(             # Add a box around the plot
  #  color = "grey", size = 0.5, fill = NA            # Black border, no fill
  #)) +
guides("none")



## ----fig.asp = 1, fig.width = 12----------------------------------------------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Define vertical reference lines for key historical periods
# ------------------------------------------------------------------------------
highlighted_dates <- c(750, 1050, 1450, 1750)

# ------------------------------------------------------------------------------
# Combine clustering outputs from all lakes into a single data frame
# ------------------------------------------------------------------------------
data_combined <- ls_all_hclusts_complete_paleo %>%
  map_dfr(as.data.frame)

# ------------------------------------------------------------------------------
# Compute time coverage per lake and order lakes by coverage length
# ------------------------------------------------------------------------------
lake_ranges <- data_combined %>%
  group_by(lake) %>%
  summarise(
    min_age = min(age_ce, na.rm = TRUE),
    max_age = max(age_ce, na.rm = TRUE),
    data_range = max_age - min_age,
    .groups = "drop"
  ) %>%
  arrange(desc(data_range))

lake_levels <- lake_ranges$lake  # Preserve lake order for plotting

# ------------------------------------------------------------------------------
# Set lake as ordered factor and cluster as categorical
# ------------------------------------------------------------------------------
data_combined <- data_combined %>%
  mutate(
    lake = factor(lake, levels = lake_levels),
    amd_clusts = as.factor(amd_clusts)
  )

# ------------------------------------------------------------------------------
# Define shaded rectangles to mask periods without data
# ------------------------------------------------------------------------------
rect_data <- lake_ranges %>%
  select(lake, min_age, max_age) %>%
  mutate(
    lake = factor(lake, levels = lake_levels),
    xmin1 = 0,
    xmax1 = min_age,
    xmin2 = max_age,
    xmax2 = 2010
  ) %>%
  pivot_longer(c(xmin1, xmin2), names_to = "rect_part", values_to = "xmin") %>%
  pivot_longer(c(xmax1, xmax2), names_to = "rect_part2", values_to = "xmax") %>%
  filter(substr(rect_part, 5, 5) == substr(rect_part2, 5, 5)) %>%
  select(lake, xmin, xmax) %>%
  mutate(ymin = -Inf, ymax = Inf)

# ------------------------------------------------------------------------------
# Keep only lakes with ≥2 trophic clusters
# ------------------------------------------------------------------------------
lake_valid_clusters <- data_combined %>%
  distinct(lake, amd_clusts) %>%
  count(lake) %>%
  filter(n >= 2) %>%
  pull(lake)

# ------------------------------------------------------------------------------
# Final filtering: remove sparse groups and ensure valid lake set
# ------------------------------------------------------------------------------
data_combined_clean <- data_combined %>%
  filter(
    lake %in% lake_valid_clusters,
    !is.na(age_ce), !is.na(amd_clusts),
    is.finite(age_ce)
  ) %>%
  group_by(lake, amd_clusts) %>%
  filter(n() >= 3, sd(age_ce, na.rm = TRUE) > 1) %>%
  ungroup()

rect_data_clean <- rect_data %>%
  filter(lake %in% lake_valid_clusters)

# ------------------------------------------------------------------------------
# Build final plot
# ------------------------------------------------------------------------------
plot_density_fclusts_time_lines_1B <- 
  ggplot(data_combined_clean, aes(x = age_ce, fill = amd_clusts)) +
  facet_wrap(~ lake, ncol = 1, scales = "fixed", strip.position = "right") +
  geom_density(adjust = 1.5, position = "fill", colour = "white") +
  geom_rect(
    data = rect_data_clean,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "white", colour = "white"
  ) +
  geom_vline(xintercept = highlighted_dates, colour = "white", linetype = "dashed") +
  scale_fill_viridis_d(
    name = "Community Trophic \nstructures", option = "viridis", direction = -1
  ) +
  lims(x = c(0, 2010)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  theme_minimal() +
  theme(
    plot.background = element_blank(),
    plot.title = element_text(size = 12),
    legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
    legend.justification = c(0, 1),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
    strip.text = element_text(size = 10),
    plot.margin = unit(c(2, 1, 1, 1), "cm")
  ) +
  labs(
    x = "Age (CE)",
    y = "Proportion of Community Trophic Structures",
    title = "Proportions of Community Trophic Structures Over Time by Lake"
  )

# ------------------------------------------------------------------------------
# Render final plot (warning suppressed due to known density dropouts)
# ------------------------------------------------------------------------------
suppressWarnings(print(plot_density_fclusts_time_lines_1B))



## ----fig.asp = 0.5, fig.width = 12--------------------------------------------------------------------------------------------------------------------------------------------------------------

# Create the plot
plot_stacked_fclusts_time_lines <- 
  ls_all_hclusts_complete_paleo %>%
  map(~.x %>% as.data.frame()) %>% 
  bind_rows() %>%
  ggplot(aes(x = age_ce, fill = amd_clusts)) +
  geom_density(adjust = 1.5, position = "fill", colour = "white") +
  scale_fill_viridis_d(name = "Community Trophic \nstructures", option = "viridis", direction = -1) +
  lims(x = c(0, 2010)) +
  # Add vertical lines for all highlighted dates
  geom_vline(data = data.frame(x = highlighted_dates), aes(xintercept = x), colour = "white") +
  # Add labels for the dates
  geom_text(
    data = data.frame(x = highlighted_dates, y = 1.05, label = as.character(highlighted_dates)),
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    colour = "black",
    size = 3,
    fontface = "bold"
  ) +
  theme(plot.background = NULL) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12),
    legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
    legend.justification = c(0, 1),
    # legend.position = c(0, 1),
    axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
    panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Age (CE)",
    y = "Proportion of Community Trophic Structures",
    title = "Combined proportions of Community Trophic Structures over Time"
  )

plot_stacked_fclusts_time_lines


## ----fig.asp = 0.5, fig.width = 12--------------------------------------------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Compute ordination using PCoA on Bray-Curtis dissimilarity
# ------------------------------------------------------------------------------
ordinations_paleo <- norm_abund %>%
  as.data.frame() %>%
  select(-c(lake:age_ce)) %>%                        # Drop metadata columns
  vegdist(method = "bray", na.rm = TRUE) %>%         # Bray-Curtis distance
  cmdscale() %>%                                     # PCoA
  scores() %>%
  as.data.frame()

# ------------------------------------------------------------------------------
# Fit environmental vectors (alternative functional groups) to ordination
# ------------------------------------------------------------------------------
fit_env_paleo <- ordinations_paleo %>%
  envfit(
    env = select(norm_abund, alternative_functional_groups_guilds_globi_merged),
    scaling = "sites"
  ) %>%
  scores("vectors") %>%
  as.data.frame() %>%
  mutate(
    Dim1 = Dim1 * 0.75,
    Dim2 = Dim2 * 0.75,
    groups = rownames(.)
  )

# ------------------------------------------------------------------------------
# Compute convex hulls (polygons) for each trophic structure (cluster)
# ------------------------------------------------------------------------------
hull_paleo <- bind_cols(
    ordinations_paleo,
    select(norm_abund, alternative_functional_groups_guilds_globi_merged),
    select(global_matrix_merged_paleo, DCA1:amd_clusts)
  ) %>%
  group_by(amd_clusts) %>%
  dplyr::slice(chull(Dim1, Dim2))  # Convex hull per cluster

# ------------------------------------------------------------------------------
# Rename cluster levels for better readability (CTS1–CTS5)
# ------------------------------------------------------------------------------
global_matrix_merged_paleo$amd_clusts <- factor(
  global_matrix_merged_paleo$amd_clusts,
  levels = 1:5,
  labels = c("CTS1", "CTS2", "CTS3", "CTS4", "CTS5")
)

hull_paleo$amd_clusts <- factor(
  hull_paleo$amd_clusts,
  levels = 1:5,
  labels = c("CTS1", "CTS2", "CTS3", "CTS4", "CTS5")
)

# ------------------------------------------------------------------------------
# Assemble full ordination plot with polygons, points, and fitted vectors
# ------------------------------------------------------------------------------
ordinations_fclusts_1 <- ordinations_paleo %>%
  bind_cols(
    select(norm_abund, -alternative_functional_groups_guilds_globi_merged),
    .,
    select(global_matrix_merged_paleo, DCA1:amd_clusts)
  ) %>%
  ggplot(aes(x = Dim1, y = Dim2)) +
  
  # Polygons for each trophic structure
  geom_polygon(
    data = hull_paleo,
    aes(fill = as.factor(amd_clusts)),
    alpha = 0.2
  ) +
  scale_fill_viridis_d(option = "viridis", direction = -1, guide = "none") +
  ggnewscale::new_scale_fill() +
  
  # Points colored by cluster
  geom_point(aes(colour = as.factor(amd_clusts)), size = 3) +
  scale_colour_viridis_d(
    option = "viridis", direction = -1,
    name = "Community Trophic Structures"
  ) +
  ggnewscale::new_scale_colour() +

  # Environmental vectors
  geom_segment(
    data = fit_env_paleo,
    aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
    arrow = arrow(length = unit(0.1, "cm"))
  ) +

  # Environmental labels
  ggrepel::geom_label_repel(
    data = fit_env_paleo,
    aes(x = Dim1, y = Dim2, label = groups),
    max.overlaps = 12,
    size = 3,
    direction = "both",
    segment.size = 0.25
  ) +

  # Styling
  theme_bw() +
  labs(title = "") +
  theme(
    plot.title = element_text(size = 12),
    legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = NULL,
    axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
    panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
    panel.grid.minor = element_blank()
  ) +
  guides(fill = "none")

# ------------------------------------------------------------------------------
# Render plot
# ------------------------------------------------------------------------------
ordinations_fclusts_1



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# ----------------------------
# Data Preparation
# ----------------------------

# Combine global matrix with normalized abundance data (excluding lake to age_ce columns)
datRF_paleo <- bind_cols(
  global_matrix_merged_paleo,
  select(norm_abund, -c(lake:age_ce))
) %>% 
  select(amd_clusts:last_col()) %>%     # Keep columns from amd_clusts to the last
  select(-total_nspp_by_lake_core)      # Remove total species richness column

# ----------------------------
# Random Forest Classification
# ----------------------------

# Convert cluster variable to factor for classification
datRF_paleo$amd_clusts <- as.factor(datRF_paleo$amd_clusts)

# Display class distribution
table(datRF_paleo$amd_clusts)

# Get the minimum class size (for potential balancing use)
min <- as.numeric(min(table(datRF_paleo$amd_clusts)))

# Fit random forest classifier with 1000 trees and variable importance enabled
RFfit_paleo <- randomForest(amd_clusts ~ ., data = datRF_paleo, ntree = 1000, importance = TRUE)

# ----------------------------
# 1.1 Out-of-Bag (OOB) Prediction Accuracy
# ----------------------------

# Predict using OOB samples (internal validation)
pred_paleo <- predict(RFfit_paleo)

# Append predictions to the original data
datRF_paleo <- cbind(pred_paleo, datRF_paleo)

# Confusion matrix: observed vs predicted clusters
tabrp_paleo <- table(datRF_paleo$amd_clusts, datRF_paleo$pred_paleo)

# Print row-wise classification accuracy (%)
round(100 * (prop.table(tabrp_paleo, 1)), 0)

# Print overall classification accuracy
paste("Classification accuracy:", round((100 * (sum(diag(tabrp_paleo)) / sum(tabrp_paleo))), 0), "%")

# Compute Cohen’s Kappa for classification agreement
kappa <- kappa2(datRF_paleo[, c(1, 2)], "equal")
round(kappa$value, 2)  # Interpret as moderate if between 0.4 and 0.6

# Diagonal: per-class correct classification %
diag(round(100 * (prop.table(tabrp_paleo, 1)), 0))

# ----------------------------
# Variable Importance and Partial Dependence Profiles (PDPs)
# ----------------------------

# Use DALEX to explain model behavior
explainer_paleo <- explain(model = RFfit_paleo, datRF_paleo[, -c(1:2)], datRF_paleo$amd_clusts)

# Compute and display variable importance via DALEX
imp_paleo <- variable_importance(explainer_paleo)
imp_paleo
plot(imp_paleo)

# Compute and plot PDPs (marginal effects of predictors)
pdp10_paleo <- model_profile(explainer_paleo)
plot(pdp10_paleo)

# ----------------------------
# Decision Tree Visualization
# ----------------------------

# Re-confirm amd_clusts is treated as a factor
datRF_paleo <- datRF_paleo %>%
  mutate(amd_clusts = as.factor(amd_clusts))

# Fit a decision tree using rpart (excluding predictions and response itself)
decision_tree_paleo <- rpart(
  datRF_paleo$amd_clusts ~ .,
  data = select(datRF_paleo, -c(pred_paleo, amd_clusts)),
  method = "class",
  cp = 0.035
)

# Display complexity parameter table
printcp(decision_tree_paleo)

# --------------------------------------------
# Save the decision tree as a high-resolution PNG image
# --------------------------------------------
png("outputs/decision_tree_paleo.png", width = 2000, height = 1400, res = 300)
rpart.plot::rpart.plot(
  decision_tree_paleo, 
  type = 5, 
  extra = 2, 
  legend.x = NA, 
  box.palette = as.list(viridis::viridis(n = ls_opt_clusts, direction = -1, alpha = 0.75))
)
dev.off()



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

lakes_col <- c(viridis::viridis(9, option = "plasma", alpha = 0.8))
lakes_col_10 <- c("black", viridis::viridis(10, option = "plasma", alpha = 0.8))

b_ls_dca_diatoms <-
  b_ls_dca_diatoms %>%
  dplyr::filter(!is.na(age_ce)) %>% 
  dplyr::filter(!lake %in% c("Fogo", "Furnas" ))

# global model diatoms
modGI_diatoms <- 
  b_ls_dca_diatoms %>%
  mutate(lake = as.factor(lake)) %>%
  gam(
    DCA1 ~
      s(age_ce, bs = "tp", k = 30) +
      s(age_ce, by = lake, k = 8, m = 1, bs = "tp") +
      s(lake, bs = "re", k = 8),
    data = ., method = "REML"
  )

# Output p.Time.GI, p.time.derivative
mod_diatoms <- modGI_diatoms

# Model diagnostics
summary(mod_diatoms)
k.check(mod_diatoms, subsample=5000, n.rep=400)
gam.check(mod_diatoms)
appraise(mod_diatoms)
gratia::draw(mod_diatoms, residuals = T)

# get the model data to plot
# evaluate the smooths
smooth_diatoms <- smooth_estimates(mod_diatoms) %>%
  add_confint()

# add partial residuals to data
df_partial_residuals_diatoms <- 
  b_ls_dca_diatoms %>%
  add_partial_residuals(., model = mod_diatoms)

# add label
sum1 <- summary(mod_diatoms)
lab1 <- paste("AdjRsq = ", round(sum1$r.sq,2),
              "\nDevExpl = ", round(sum1$dev.expl,2))

plot_time_GI_diatoms <- 
  smooth_diatoms %>%
  filter(.smooth == "s(age_ce)") %>%
  ggplot() +
  geom_rug(aes(x = age_ce),
    data = df_partial_residuals_diatoms,
    sides = "b", length = grid::unit(0.02, "npc")
  ) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, x = age_ce),
    alpha = 0.2
  ) +
  # geom_point(aes(x = age_ce, y = `s(age_ce)`, color = lake),
  #           data = df_partial_residuals_diatoms, alpha = 0.5, size = 1) +
  geom_line(aes(x = age_ce, y = .estimate), lwd = 1) +
  labs(y = "DCA1", title = "s(age_ce) GI model") +
  scale_color_manual("Time-Interval") +
  theme_classic(base_size = 12) +
  scale_x_continuous(n.breaks = 5, limits = c(0, 2020)) +
  lims(y = c(-2, 2)) +
  labs(title = NULL) +
  theme(panel.grid.minor.y = element_blank()) +
  theme(legend.position = "none") +
    annotate("text_npc", npcx = 0.05, npcy = .95, label = lab1, cex = 3.5)

plot_time_GI_diatoms



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

dev_mod_diatoms <- derivatives(mod_diatoms, term = "s(age_ce)")

plot_time_dev_diatoms <-
  dev_mod_diatoms %>%
  ggplot() +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, x = age_ce),
    alpha = 0.2
  ) +
  geom_line(aes(x = age_ce, y = .derivative, color = (.lower_ci) > 0 | (.upper_ci) < 0, group = 1),
    lwd = 1, show.legend = FALSE
  ) +
  scale_x_continuous(n.breaks = 5, limits = c(0, 2020)) +
  scale_y_continuous(name = "derivative", n.breaks = 3, limits = c(-0.010, 0.010)) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic(base_size = 12) +
  lims(x = c(0, 2020)) +
  lims(y = c(-0.025, 0.025)) +
  theme(panel.grid.minor.y = element_blank())

plot_time_dev_diatoms



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

b_ls_dca_chiro <-
  b_ls_dca_chiro %>%
  dplyr::filter(!is.na(age_ce)) %>% 
  dplyr::filter(!lake %in% c("Fogo", "Furnas" ))

# global model chiro
modGI_chiro <- 
  b_ls_dca_chiro %>%
  mutate(lake = as.factor(lake)) %>%
  gam(
    DCA1 ~
      s(age_ce, bs = "tp", k = 30) +
      s(age_ce, by = lake, k = 8, m = 1, bs = "tp") +
      s(lake, bs = "re", k = 8),
    data = ., method = "REML"
  )

# Output p.Time.GI, p.time.derivative
mod_chiro <- modGI_chiro

# Model diagnostics
summary(mod_chiro)
k.check(mod_chiro, subsample=5000, n.rep=400)
gam.check(mod_chiro)
appraise(mod_chiro)
gratia::draw(mod_chiro, residuals = T)

# get the model data to plot
# evaluate the smooths
smooth_chiro <- smooth_estimates(mod_chiro) %>%
  add_confint()

# add partial residuals to data
df_partial_residuals_chiro <- 
  b_ls_dca_chiro %>%
  add_partial_residuals(., model = mod_chiro)

# add label
sum1 <- summary(mod_chiro)
lab1 <- paste("AdjRsq = ", round(sum1$r.sq,2),
              "\nDevExpl = ", round(sum1$dev.expl,2))

plot_time_GI_chiro <- 
  smooth_chiro %>%
  filter(.smooth == "s(age_ce)") %>%
  ggplot() +
  geom_rug(aes(x = age_ce),
    data = df_partial_residuals_chiro,
    sides = "b", length = grid::unit(0.02, "npc")
  ) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, x = age_ce),
    alpha = 0.2
  ) +
  # geom_point(aes(x = age_ce, y = `s(age_ce)`, color = lake),
  #           data = df_partial_residuals_chiro, alpha = 0.5, size = 1) +
  geom_line(aes(x = age_ce, y = .estimate), lwd = 1) +
  labs(y = "DCA1", title = "s(age_ce) GI model") +
  scale_color_manual("Time-Interval", values = lakes_col) +
  theme_classic(base_size = 12) +
  scale_x_continuous(n.breaks = 5, limits = c(0, 2020)) +
  lims(y = c(-2, 2)) +
  labs(title = NULL) +
  theme(panel.grid.minor.y = element_blank()) +
  theme(legend.position = "none") +
    annotate("text_npc", npcx = 0.05, npcy = .95, label = lab1, cex = 3.5)

plot_time_GI_chiro



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

dev_mod_chiro <- derivatives(mod_chiro, term = "s(age_ce)")

plot_time_dev_chiro <-
  dev_mod_chiro %>%
  ggplot() +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, x = age_ce),
    alpha = 0.2
  ) +
  geom_line(aes(x = age_ce, y = .derivative, color = (.lower_ci) > 0 | (.upper_ci) < 0, group = 1),
    lwd = 1, show.legend = FALSE
  ) +
  scale_x_continuous(n.breaks = 5, limits = c(0, 2020)) +
  scale_y_continuous(name = "derivative", n.breaks = 3, limits = c(-0.010, 0.010)) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic(base_size = 12) +
  lims(x = c(0, 2020)) +
  lims(y = c(-0.025, 0.025)) +
  theme(panel.grid.minor.y = element_blank())

plot_time_dev_chiro


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

b_ls_dca_fgroup <-
  b_ls_dca_fgroup %>%
  dplyr::filter(!is.na(age_ce)) %>% 
  dplyr::filter(!lake %in% c("Fogo", "Furnas" ))

# global model fgroup
modGI_fgroup <- 
  b_ls_dca_fgroup %>%
  mutate(lake = as.factor(lake)) %>%
  gam(
    DCA1 ~
      s(age_ce, bs = "tp", k = 30) +
      s(age_ce, by = lake, k = 8, m = 1, bs = "tp") +
      s(lake, bs = "re", k = 8),
    data = ., method = "REML"
  )

# Output p.Time.GI, p.time.derivative
mod_fgroup <- modGI_fgroup

# Model diagnostics
summary(mod_fgroup)
k.check(mod_fgroup, subsample=5000, n.rep=400)
gam.check(mod_fgroup)
appraise(mod_fgroup)
gratia::draw(mod_fgroup, residuals = T)

# get the model data to plot
# evaluate the smooths
smooth_fgroup <- smooth_estimates(mod_fgroup) %>%
  add_confint()

# add partial residuals to data
df_partial_residuals_fgroup <- 
  b_ls_dca_fgroup %>%
  add_partial_residuals(., model = mod_fgroup)

# add label
sum1 <- summary(mod_fgroup)
lab1 <- paste("AdjRsq = ", round(sum1$r.sq,2),
              "\nDevExpl = ", round(sum1$dev.expl,2))

plot_time_GI_fgroup <- 
  smooth_fgroup %>%
  filter(.smooth == "s(age_ce)") %>%
  ggplot() +
  geom_rug(aes(x = age_ce),
    data = df_partial_residuals_fgroup,
    sides = "b", length = grid::unit(0.02, "npc")
  ) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, x = age_ce),
    alpha = 0.2
  ) +
  # geom_point(aes(x = age_ce, y = `s(age_ce)`, color = lake),
  #           data = df_partial_residuals_fgroup, alpha = 0.5, size = 1) +
  geom_line(aes(x = age_ce, y = .estimate), lwd = 1) +
  labs(y = "DCA1", title = "s(age_ce) GI model") +
  scale_color_manual("Time-Interval", values = lakes_col) +
  theme_classic(base_size = 12) +
  scale_x_continuous(n.breaks = 5, limits = c(0, 2020)) +
  lims(y = c(-2, 2)) +
  labs(title = NULL) +
  theme(panel.grid.minor.y = element_blank()) +
  theme(legend.position = "none") +
    annotate("text_npc", npcx = 0.05, npcy = .95, label = lab1, cex = 3.5)

plot_time_GI_fgroup



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

dev_mod_fgroup <- derivatives(mod_fgroup, term = "s(age_ce)")

plot_time_dev_fgroup <-
  dev_mod_fgroup %>%
  ggplot() +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, x = age_ce),
    alpha = 0.2
  ) +
  geom_line(aes(x = age_ce, y = .derivative, color = (.lower_ci) > 0 | (.upper_ci) < 0, group = 1),
    lwd = 1, show.legend = FALSE
  ) +
  scale_x_continuous(n.breaks = 5, limits = c(0, 2020)) +
  scale_y_continuous(name = "derivative", n.breaks = 3, limits = c(-0.010, 0.010)) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic(base_size = 12) +
  lims(x = c(0, 2020)) +
  lims(y = c(-0.025, 0.025)) +
  theme(panel.grid.minor.y = element_blank())

plot_time_dev_fgroup


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Rescaling Environmental variables
# ------------------------------------------------------------------------------

df_env <- 
  ls_df_env_wide_full %>% 
  bind_rows() %>% 
  dplyr::filter(!is.na(age_ce)) %>% 
  dplyr::filter(!lake %in% c("Fogo", "Furnas" )) %>% 
  group_by(lake) %>%
  mutate(across(c(tc, tn, toc_tn, d13c, d15n), ~ scale(.) %>% as.numeric())) %>%
  ungroup()



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------
# Load and clean pollen data
# ----------------------------------------------------------

pollen <- pollen %>% 
  rename(age_bp = age) %>% 
  mutate(age_ce = 1950 - age_bp) %>% 
  mutate(ID = paste0(sitename, "_", age_ce))

df_1 <- pollen

# Filter out aquatic and unidentified taxa, then group by ecological group
df_wk <- df_1 %>%
  filter(!ecologicalgroup %in% c("UNID", "AQVP", "AQBR")) %>%
  group_by(ID, sitename, age_ce, group) %>%
  summarise(value = sum(value), .groups = 'drop') %>%
  group_by(ID) %>%
  mutate(percentage = value / sum(value) * 100) %>%
  ungroup()

# Calculate group-wise mean percentages and anomalies
df_wk <- df_wk %>%
  group_by(sitename, group) %>%
  mutate(mean_percentage = mean(percentage)) %>%
  ungroup() %>%
  mutate(percent_anomaly = percentage - mean_percentage) %>%
  select(-mean_percentage) %>%
  droplevels()

sites_col <- scico(n_distinct(df_wk$sitename), begin = 0, end = 0.95, palette = 'roma')
gr_col <- scico(n_distinct(df_wk$group), begin = 0, end = 0.95, palette = 'nuuk')

# ----------------------------------------------------------
# DCA Analysis
# ----------------------------------------------------------

exclude_sites <- c("Lomba", "Rasa")
df_wk_dca <- df_1 %>%
  filter(!ecologicalgroup %in% c("UNID"), age_ce > 0, !sitename %in% exclude_sites) %>%
  droplevels()

results_list <- lapply(sort(unique(df_wk_dca$sitename)), function(site) {
  df.sp <- df_wk_dca %>%
    filter(sitename == site) %>%
    group_by(ID, taxa) %>%
    summarise(value = sum(value), .groups = 'drop') %>%
    pivot_wider(names_from = taxa, values_from = value, values_fill = 0) %>%
    column_to_rownames("ID") %>%
    select(where(~ any(. != 0)))
  
  df.spe <- decostand(df.sp, method = "hellinger")
  df.DCA <- decorana(df.spe) %>%
    scores(choices = 1) %>%
    as.data.frame() %>%
    rownames_to_column("ID")
  
  return(df.DCA)
})

DCA.res <- Reduce(full_join, results_list) %>%
  separate(ID, into = c("sitename", "age_ce"), sep = "_", convert = TRUE) %>%
  arrange(sitename, age_ce) %>%
  mutate(DCA1 = case_when(
    sitename == "Azul" ~ -DCA1,
    sitename == "Gingal" ~ -DCA1,
    TRUE ~ DCA1
  ))

DF_to <- DCA.res %>% mutate(sitename = as.factor(sitename))

# ----------------------------------------------------------
# HGAM on DCA
# ----------------------------------------------------------

library(mgcv)
mod_to <- gam(DCA1 ~ s(age_ce, bs = "tp", k = 25) +
                s(age_ce, by = sitename, bs = "tp", k = 15, m = 1) +
                s(sitename, bs = "re"),
              data = DF_to, method = "REML")

summary(mod_to)
k.check(mod_to, n.rep = 500)
appraise(mod_to)

# Plot smoothed global trend and derivative
sm_to <- smooth_estimates(mod_to) %>% add_confint()
df <- DF_to %>% add_partial_residuals(mod_to)
lakes_col <- scico(length(levels(DF_to$sitename)), begin = 0, end = 0.95, palette = 'roma')

p.DCA.gl <- sm_to %>%
  filter(.smooth == "s(age_ce)") %>%
  ggplot() +
  geom_rug(aes(x = age_ce), data = df, sides = "b", length = unit(0.02, "npc")) +
  geom_ribbon(aes(x = age_ce, ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2) +
  geom_point(aes(x = age_ce, y = `s(age_ce)`, color = sitename), data = df, alpha = 0.8, size = 1.5) +
  geom_line(aes(x = age_ce, y = .estimate), linewidth = 1) +
  labs(y = "DCA axis-1", x = "Age (CE)") +
  scale_color_manual("Lake", values = lakes_col) +
  theme_minimal(base_size = 9) +
  scale_x_continuous(n.breaks = 10) +
  theme(panel.grid.minor.y = element_blank(), legend.position = "bottom")

mod_to_dev <- derivatives(mod_to, select = "s(age_ce)", eps = 1, level = 0.9, interval = "simultaneous", seed = 22)
color_mapping <- scico(2, begin = 0, end = 0.9, palette = 'vik')

# ----------------------------------------------------------
# HGAM on representative groups
# ----------------------------------------------------------

pollen_clean <- pollen %>%
  filter(elementtype == "pollen", !is.na(value)) %>%
  filter(!sitename %in% c("Lomba", "Rasa"))

sample_totals <- pollen_clean %>%
  group_by(sitename, age_ce) %>%
  summarise(sample_total = sum(value, na.rm = TRUE), .groups = "drop")

lake_avg <- sample_totals %>%
  group_by(sitename) %>%
  summarise(avg_lake_total = mean(sample_total, na.rm = TRUE), .groups = "drop")

pollen_normalized <- pollen_clean %>%
  left_join(sample_totals, by = c("sitename", "age_ce")) %>%
  left_join(lake_avg, by = "sitename") %>%
  mutate(
    scaling_factor = avg_lake_total / sample_total,
    value_scaled = value * scaling_factor
  )

# --- Calculate percent composition by group ---

pollen_percent_group <- pollen_normalized %>%
  group_by(sitename, age_ce, group) %>%
  summarise(value_scaled = sum(value_scaled, na.rm = TRUE), .groups = "drop") %>%
  group_by(sitename, age_ce) %>%
  mutate(
    total_scaled = sum(value_scaled),
    percent = (value_scaled / total_scaled) * 100
  ) %>%
  filter(!is.na(percent)) %>%
  mutate(group = factor(group), sitename = factor(sitename)) %>%
  ungroup()

pollen_percent_group <- pollen_percent_group %>%
  mutate(sitename = factor(sitename))  # ensure factor type

# --- GAM fitting function ---

fit_gam_group <- function(df) {
  if (!"sitename" %in% names(df)) {
    stop("Missing `sitename` column in group: ", unique(df$group))
  }
  
  if (n_distinct(df$sitename) < 2) {
    warning("Only one site found in group: ", unique(df$group), ". Fitting a basic model.")
    return(
      gam(percent ~ s(age_ce, k = 10), data = df, method = "REML")
    )
  }

  gam(
    percent ~ 
      s(age_ce, bs = "tp", k = 30) +
      s(age_ce, by = sitename, bs = "tp", k = 9, m = 1) +
      s(sitename, bs = "re"),
    data = df,
    method = "REML"
  )
}

# --- Fit GAM for each pollen group ---

gam_models <- pollen_percent_group %>%
  group_split(group) %>%
  set_names(map_chr(., ~ as.character(unique(.x$group)))) %>%
  map(fit_gam_group)

# --- Plotting function per group ---

plot_gam_group <- function(group_name, mod) {
  df_aug <- augment(mod, type.predict = "response")
  sum1 <- summary(mod)
  lab1 <- paste("AdjRsq = ", round(sum1$r.sq, 2),
                "\nDevExpl = ", round(sum1$dev.expl, 2))
  
  ggplot(df_aug, aes(x = age_ce)) +
    geom_point(aes(y = percent), color = "black", alpha = 0.3, size = 1) +
    geom_smooth(aes(y = .fitted), method = "loess", se = TRUE, color = "steelblue", linewidth = 1) +
    labs(title = group_name, y = "% Abundance") +
    scale_x_continuous(limits = c(0, 2020), n.breaks = 5) +
    annotate("text_npc", npcx = 0.05, npcy = 0.95, label = lab1, size = 3.5, hjust = 0, vjust = 1) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
}

# --- Generate plots ---

plots_all <- map2(names(gam_models), gam_models, plot_gam_group)

final_plot <- wrap_plots(plots_all, ncol = 2) +
  plot_annotation(
    title = "% of Pollen Groups Over Time (GAM fit)",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

# --- Display plot ---

final_plot

# --- Compare tree vs herb pollen groups ---

tree_data <- pollen_percent_group %>% filter(group == "Tree")
herb_data <- pollen_percent_group %>% filter(group == "Herb")

mod_tree <- fit_gam_group(tree_data)
mod_herb <- fit_gam_group(herb_data)

aug_tree <- augment(mod_tree, type.predict = "response") %>%
  mutate(group = "Tree")

aug_herb <- augment(mod_herb, type.predict = "response") %>%
  mutate(group = "Herb")

combined_df_pollen <- bind_rows(aug_tree, aug_herb)



## ----fig.width = 12, fig.asp = 1.2--------------------------------------------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Interpolate vegetation smooth (HGAM output)
# ------------------------------------------------------------------------------

interp_years <- 0:2010

# Extract the global smooth only (sitename == NA)
sm_interp <- sm_to %>%
  filter(is.na(sitename)) %>%
  select(age_ce, estimate = .estimate, lower_ci = .lower_ci, upper_ci = .upper_ci)

# Interpolate all three columns at once using purrr and bind together
interpolated_df <- interp_years %>%
  tibble(age_ce = .) %>%
  mutate(
    estimate = approx(x = sm_interp$age_ce, y = sm_interp$estimate, xout = age_ce, method = "linear")$y,
    lower_ci = approx(x = sm_interp$age_ce, y = sm_interp$lower_ci, xout = age_ce, method = "linear")$y,
    upper_ci = approx(x = sm_interp$age_ce, y = sm_interp$upper_ci, xout = age_ce, method = "linear")$y
  )

# ------------------------------------------------------------------------------
# Prepare and interpolate NAO Hernandez
# ------------------------------------------------------------------------------

colnames(nao_hernandez)[2] <- "age_ce"
colnames(nao_hernandez)[3] <- "nao_lower_ci"
colnames(nao_hernandez)[5] <- "NAO_Median_Value"
colnames(nao_hernandez)[7] <- "nao_upper_ci"

nao_hernandez <- nao_hernandez %>%
  filter(age_ce >= 0) %>%
  mutate(sign = ifelse(NAO_Median_Value >= 0, "positive", "negative"))

nao_hernandez_interpolated <- approx(
  x = nao_hernandez$age_ce,
  y = nao_hernandez$NAO_Median_Value,
  xout = interp_years,
  method = "linear"
)

nao_hernandez_interpolated <- as.data.frame(nao_hernandez_interpolated)
colnames(nao_hernandez_interpolated) <- c("age_ce", "NAO_Median_Value")

nao_binned <- nao_hernandez_interpolated %>%
  mutate(bin_50yr = floor(age_ce / 30) * 30) %>%
  group_by(bin_50yr) %>%
  summarise(avg_value = mean(NAO_Median_Value, na.rm = TRUE), .groups = "drop")

nao_hernandez <- nao_hernandez %>% mutate(age_ce = round(age_ce))

# ------------------------------------------------------------------------------
# Prepare and interpolate NH summer temperatures (Büntgen et al.)
# ------------------------------------------------------------------------------

nhst_buntgen_interpolated <- approx(
  x = nhst_buntgen$age_ce,
  y = nhst_buntgen$Rmean,
  xout = interp_years,
  method = "linear"
)

nhst_buntgen <- as.data.frame(nhst_buntgen_interpolated)
colnames(nhst_buntgen) <- c("age_ce", "Rmean")

nhst_buntgen <- nhst_buntgen %>%
  mutate(sign = ifelse(Rmean >= 0, "positive", "negative"))

nhst_buntgen_binned <- nhst_buntgen %>%
  mutate(bin_50yr = floor(age_ce / 50) * 50) %>%
  group_by(bin_50yr) %>%
  summarise(avg_value = mean(Rmean, na.rm = TRUE), .groups = "drop")

# ------------------------------------------------------------------------------
# Prepare normalized guild abundances and community matrix
# ------------------------------------------------------------------------------

norm_abund_guilds <- norm_abund_guilds %>%
  mutate(age_ce = round(age_ce))

comm_matrix <- norm_abund_guilds %>%
  select(high_profile:predator)

norm_abund_guilds <- norm_abund_guilds %>%
  select(-high_profile:-predator) %>%
  bind_cols(comm_matrix)

# ------------------------------------------------------------------------------
# Estimate temporal density and relative proportions of CTS clusters
# ------------------------------------------------------------------------------

combined_data <- ls_all_hclusts_complete_paleo %>%
  map(~as.data.frame(.x)) %>%
  bind_rows()

densities <- combined_data %>%
  group_by(amd_clusts) %>%
  do({
    dens <- density(.$age_ce, adjust = 1.5, from = 0, to = 2010, n = 5001)
    data.frame(age_ce = dens$x, density = dens$y)
  }) %>%
  ungroup()

proportions_df <- densities %>%
  group_by(age_ce) %>%
  mutate(
    total_density = sum(density),
    proportion = density / total_density
  ) %>%
  ungroup()

proportions_long <- proportions_df %>%
  mutate(age_ce = round(age_ce)) %>%
  filter(age_ce >= 0, age_ce <= 2010) %>%
  group_by(age_ce, amd_clusts) %>%
  summarise(
    average_proportion = mean(proportion, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# Convert CTS proportions to wide matrix format
# ------------------------------------------------------------------------------

prop_comm_matrix <- proportions_long %>%
  tidyr::pivot_wider(
    names_from = amd_clusts,
    values_from = average_proportion
  )

# ------------------------------------------------------------------------------
# Prepare and align datasets using 30-year windows
# ------------------------------------------------------------------------------

joined_df <- cbind(age_ce = norm_abund_guilds$age_ce, comm_matrix) %>%
  left_join(nao_hernandez_interpolated, by = "age_ce") %>%
  left_join(interpolated_df, by = "age_ce") %>%
  left_join(nhst_buntgen, by = "age_ce") %>%
  left_join(prop_comm_matrix, by = "age_ce")

joined_df_30yr <- joined_df %>%
  filter(age_ce > 0) %>%
  mutate(age_bin = floor(age_ce / 30) * 30) %>%  # define 30-year bins
  group_by(age_bin) %>%
  summarise(
    across(where(is.numeric) & !all_of("age_ce"), mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(age_ce = age_bin)  # unify time axis for modeling

# ------------------------------------------------------------------------------
# Scale environmental variables
# ------------------------------------------------------------------------------

joined_df_30yr <- joined_df_30yr %>%
  mutate(across(NAO_Median_Value:Rmean, scale))



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Prepare long-format data for GAM fitting
# ------------------------------------------------------------------------------

g_norm_abund_guilds <- 
  norm_abund_guilds %>% 
  pivot_longer(
    cols = `high_profile`:predator,
    names_to = "guild",
    values_to = "rel_abund"
  )

# Define list of trophic guilds
guilds <- c("high_profile", "low_profile", "motile", "euplanktonic",
            "algivore", "detritivore", "plantivore", "predator")

# ------------------------------------------------------------------------------
# Fit a GAM model to each trophic guild
# ------------------------------------------------------------------------------

gam_models <- guilds %>%
  set_names() %>%
  map(function(g) {
    df <- g_norm_abund_guilds %>%
      filter(guild == g) %>%
      mutate(lake = as.factor(lake))
    
    gam(
      rel_abund ~
        s(age_ce, bs = "tp", k = 30) +                                # Global smooth
        s(age_ce, by = lake, bs = "tp", k = 14, m = 1) +              # Lake-specific smooths
        s(lake, k = 14, bs = "re"),                                   # Random intercepts for lake
      data = df,
      method = "REML"
    )
  })

# ------------------------------------------------------------------------------
# Run gam.check on each model to assess k-index adequacy
# ------------------------------------------------------------------------------

gam_models_checked <- map(gam_models, mgcv::gam.check)

# ------------------------------------------------------------------------------
# Extract k-check tables into a tidy tibble
# ------------------------------------------------------------------------------

extract_kcheck_table <- function(model) {
  txt <- capture.output(gam.check(model, rep = 1000, verbose = FALSE))
  
  start_line <- grep("Basis dimension \\(k\\) checking results", txt)
  if (length(start_line) == 0) return(NULL)

  table_lines <- txt[(start_line + 2):length(txt)]
  table_lines <- table_lines[table_lines != "" & !grepl("^---", table_lines)]

  # Parse each line to extract values
  kcheck_df <- map_dfr(table_lines, function(line) {
    parts <- strsplit(trimws(line), "\\s+")[[1]]
    n_parts <- length(parts)
    smooth_term <- paste(parts[1:(n_parts - 4)], collapse = " ")
    
    tibble(
      smooth_term = smooth_term,
      k_prime     = as.numeric(parts[n_parts - 3]),
      edf         = as.numeric(parts[n_parts - 2]),
      k_index     = as.numeric(parts[n_parts - 1]),
      p_value     = suppressWarnings(as.numeric(parts[n_parts]))
    )
  })

  return(kcheck_df)
}

# ------------------------------------------------------------------------------
# Summarize all k-check results into one tibble
# ------------------------------------------------------------------------------

gam_kcheck_summary <- imap_dfr(gam_models, function(model, name) {
  kcheck_df <- extract_kcheck_table(model)
  
  if (is.null(kcheck_df)) {
    return(tibble(
      model_name  = name,
      smooth_term = NA,
      k_prime     = NA,
      edf         = NA,
      k_index     = NA,
      p_value     = NA
    ))
  }
  
  kcheck_df %>%
    mutate(model_name = name, .before = 1)
})



## ----fig.width = 12, fig.asp = 1----------------------------------------------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Prepare community and environmental matrices
# ------------------------------------------------------------------------------

vars_to_test <- c("NAO_Median_Value", "Rmean")

comm <- joined_df_30yr %>%
  select(`high_profile`:predator)

env <- joined_df_30yr %>%
  select(c(age_ce, vars_to_test, estimate))

# ------------------------------------------------------------------------------
# Check for collinearity using VIF
# ------------------------------------------------------------------------------

check_collinearity <- function(vars, threshold = 10) {
  dummy_response <- matrix(rnorm(nrow(vars)), ncol = 1)
  rda_model <- rda(dummy_response ~ ., data = vars)
  vif_values <- vif.cca(rda_model)
  cat("\nCollinearity check (VIF):\n")
  print(round(vif_values, 2))
  if (any(vif_values > threshold)) {
    warning("High collinearity detected (VIF > ", threshold, "): ",
            paste(names(vif_values[vif_values > threshold]), collapse = ", "))
  }
  return(vif_values)
}

cond_vars <- env %>% select(age_ce)
explanatory_vars <- env %>% select(-age_ce)

all_data <- bind_cols(comm, explanatory_vars, cond_vars)
complete_data <- all_data %>% filter(complete.cases(.))

comm_clean <- complete_data %>% select(colnames(comm))
explanatory_clean <- complete_data %>% select(colnames(explanatory_vars))
cond_clean <- complete_data %>% select(colnames(cond_vars))

# ------------------------------------------------------------------------------
# Fit initial CCA model and filter collinear variables
# ------------------------------------------------------------------------------

rda_full <- capscale(comm_clean ~ ., data = explanatory_clean, conditioning = cond_clean, permutations = 9999)

vif_cutoff <- 10
keep_vars <- names(explanatory_clean)

repeat {
  model <- capscale(comm_clean ~ ., data = explanatory_clean[, keep_vars, drop = FALSE])
  vif_vals <- vif.cca(model)
  if (all(vif_vals < vif_cutoff)) break
  var_to_remove <- names(vif_vals)[which.max(vif_vals)]
  cat("Removing collinear variable due to high VIF:", var_to_remove, "\n")
  keep_vars <- setdiff(keep_vars, var_to_remove)
}

rda_clean <- capscale(comm_clean ~ ., data = explanatory_clean[, keep_vars, drop = FALSE])

# ------------------------------------------------------------------------------
# Permutation tests for significance
# ------------------------------------------------------------------------------

cat("\nGlobal model significance:\n")
print(anova(rda_clean, permutations = 999))

cat("\nAxis-wise significance:\n")
print(anova(rda_clean, by = "axis", permutations = 999))

cat("\nTerm-wise (variable) significance:\n")
term_test <- anova(rda_clean, by = "term", permutations = 999)
print(term_test)

sig_terms <- rownames(term_test)[which(term_test$`Pr(>F)` < 0.05)]
sig_terms <- gsub(".*\\)", "", sig_terms)

cat("\nSignificant predictors (p < 0.05):\n")
print(sig_terms)

if (length(sig_terms) > 0) {
  cca_sig <- cca(comm_clean ~ ., data = explanatory_clean[, sig_terms, drop = FALSE])
  cat("\nRefit CCA with significant variables only:\n")
  print(summary(cca_sig))
} else {
  cat("\nNo significant predictors found at p < 0.05\n")
}

# ------------------------------------------------------------------------------
# Forward selection with ordistep
# ------------------------------------------------------------------------------

cat("\nForward selection with ordistep:\n")
null <- capscale(comm_clean ~ 1, data = explanatory_clean[, keep_vars, drop = FALSE])
full <- capscale(comm_clean ~ ., data = explanatory_clean[, keep_vars, drop = FALSE])

fwd <- ordistep(null, scope = formula(full), direction = "forward", permutations = 9999, trace = TRUE)

cat("Selected model from forward selection:\n")
print(fwd$call)

selected_formula <- formula(fwd)
selected_variables <- all.vars(selected_formula[[3]])
fwd_model <- eval(fwd$call)

biplot_scores <- scores(fwd_model, display = "bp", scaling = 2)

cat("\nCorrelations of variables with each RDA axis (biplot scores):\n")
print(round(biplot_scores, 3))

biplot_df <- as_tibble(biplot_scores, rownames = "Variable")
print(biplot_df)

# ------------------------------------------------------------------------------
# Plot RDA variable correlation heatmap
# ------------------------------------------------------------------------------

term_test <- anova(fwd_model, by = "term", permutations = 999)
term_pvals <- as.data.frame(term_test)
term_pvals$Variable <- rownames(term_pvals)
term_pvals$Variable <- gsub("`", "", term_pvals$Variable)

biplot_df <- as.data.frame(biplot_scores) %>%
  rownames_to_column("Variable") %>%
  pivot_longer(-Variable, names_to = "Axis", values_to = "Correlation") %>%
  left_join(term_pvals, by = "Variable") %>%
  mutate(Significance = ifelse(`Pr(>F)` < 0.05, "*", ""))

ggplot(biplot_df, aes(x = Axis, y = Variable, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Significance), color = "black", size = 6, vjust = 0.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Correlation") +
  labs(title = "Correlation of Variables with RDA Axes", subtitle = "* = p < 0.05 (from term-wise test)",
       x = "RDA Axis", y = "Environmental Variable") +
  theme_minimal(base_size = 13)

# ------------------------------------------------------------------------------
# Final RDA biplot figure
# ------------------------------------------------------------------------------

# Extract eigenvalues and compute explained variance
eig_vals <- fwd_model$CCA$eig
var_exp <- eig_vals / sum(eig_vals) * 100
axis1_lab <- paste0("CAP1 (", round(var_exp[1], 1), "% explained variance)")
axis2_lab <- paste0("CAP2 (", round(var_exp[2], 1), "% explained variance)")

# Get site scores and name columns
site_scores <- scores(fwd_model, display = "sites", scaling = 2)
site_df <- as.data.frame(site_scores)
colnames(site_df)[1:2] <- c("CAP1", "CAP2")  # Ensure correct column names
site_df$age_ce <- cond_clean$age_ce

# Get species scores and rename columns
species_df <- as.data.frame(scores(fwd_model, display = "species", scaling = 2))
colnames(species_df)[1:2] <- c("CAP1", "CAP2")  # Ensure correct column names
species_df <- species_df %>%
  tibble::rownames_to_column("Guild") %>%
  mutate(Label = recode(Guild,
    "high_profile" = "High Profile", "low_profile" = "Low Profile",
    "motile" = "Motile", "euplanktonic" = "Euplanktonic",
    "algivore" = "Algivore", "detritivore" = "Detritivore",
    "plantivore" = "Plantivore", "predator" = "Predator"))

# Get environmental (biplot) scores and rename columns
bp_df <- as.data.frame(scores(fwd_model, display = "bp", scaling = 2))
colnames(bp_df)[1:2] <- c("CAP1", "CAP2")  # Ensure correct column names
bp_df <- bp_df %>%
  tibble::rownames_to_column("Variable") %>%
  mutate(Label = recode(Variable,
    "estimate" = "Veg. Change",
    "NAO_Median_Value" = "NAO",
    "Rmean" = "NHST"))

# Plot
library(ggrepel)
global_rda <- ggplot() +
  geom_point(data = site_df, aes(x = CAP1, y = CAP2, color = age_ce), size = 4, alpha = 0.9) +
  geom_text_repel(data = species_df, aes(x = CAP1, y = CAP2, label = Label),
                  color = "red", size = 3.5, fontface = "italic", max.overlaps = Inf) +
  geom_segment(data = bp_df,
               aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", linewidth = 0.7) +
  geom_text_repel(data = bp_df,
                  aes(x = CAP1, y = CAP2, label = Label),
                  size = 3.5, color = "black", max.overlaps = Inf) +
  scale_color_viridis_c(option = "D", name = "Age (CE)") +
  labs(x = axis1_lab, y = axis2_lab, title = "") +
  coord_fixed() +
  theme_minimal(base_size = 14) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

# Print plot
global_rda



## ----fig.width=10, fig.asp=0.5------------------------------------------------------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------
# Phase varpart - Historical phases
# ----------------------------------------------------------

# --- Define predictor variables ---
selected_variables <- c("NAO_Median_Value", "estimate")

# --- Define phase subsets ---
phase_list <- list(
  "Phase 1" = joined_df_30yr %>% filter(age_ce < 750),
  "Phase 2" = joined_df_30yr %>% filter(age_ce >= 750 & age_ce < 1050),
  "Phase 3" = joined_df_30yr %>% filter(age_ce >= 1050 & age_ce < 1450),
  "Phase 4" = joined_df_30yr %>% filter(age_ce >= 1450 & age_ce < 1750),
  "Phase 5" = joined_df_30yr %>% filter(age_ce >= 1750)
)

# --- Function to compute pure and shared effects ---
calculate_pure_shared <- function(data, phase_name) {
  comm <- data %>% select(high_profile:predator) %>% filter(complete.cases(.))
  env  <- data %>% select(all_of(selected_variables)) %>% filter(complete.cases(.))
  n <- min(nrow(comm), nrow(env))
  comm <- comm[1:n, ]
  env  <- env[1:n, ]
  
  pure_effects <- map_dfr(selected_variables, function(var) {
    others <- setdiff(selected_variables, var)
    formula <- as.formula(
      paste0("comm ~ ", var, " + Condition(", paste(others, collapse = " + "), ")")
    )
    mod <- rda(formula = formula, data = env)
    tibble(
      phase = phase_name,
      component = ifelse(var == "estimate", "Pure Vegetation", "Pure Climate"),
      value = max(0, RsquareAdj(mod)$adj.r.squared)
    )
  })
  
  mod_full <- rda(comm ~ ., data = env)
  total_r2 <- max(0, RsquareAdj(mod_full)$adj.r.squared)
  shared_val <- total_r2 - sum(pure_effects$value)
  
  shared <- tibble(
    phase = phase_name,
    component = "Shared",
    value = max(0, shared_val)
  )
  
  bind_rows(pure_effects, shared)
}

# --- Run across all phases ---
r2_partitioned <- imap_dfr(phase_list, calculate_pure_shared)

# --- Add readable and ordered phase labels ---
r2_partitioned <- r2_partitioned %>%
  mutate(
    phase_label = factor(case_when(
      phase == "Phase 1" ~ "Phase 1\n(<750 CE)",
      phase == "Phase 2" ~ "Phase 2\n(750–1050 CE)",
      phase == "Phase 3" ~ "Phase 3\n(1050–1450 CE)",
      phase == "Phase 4" ~ "Phase 4\n(1450–1750 CE)",
      phase == "Phase 5" ~ "Phase 5\n(>1750 CE)"
    ), levels = c(
      "Phase 1\n(<750 CE)",
      "Phase 2\n(750–1050 CE)",
      "Phase 3\n(1050–1450 CE)",
      "Phase 4\n(1450–1750 CE)",
      "Phase 5\n(>1750 CE)"
    ))
  )

# --- Define muted color palette for components ---
component_colors <- c(
  "Pure Vegetation" = "#97D8C4",
  "Pure Climate" = "#F4B942",
  "Shared" = "#4059AD"
)

# --- Plot variance partitioning per phase ---
original_varpart <- ggplot(r2_partitioned, aes(x = phase_label, y = value, fill = component)) +
  geom_col(width = 0.6) +
  geom_text(
    data = r2_partitioned %>% filter(value > 0),
    aes(label = paste0(round(value * 100, 1), "%")),
    position = position_stack(vjust = 0.5),
    color = "black", size = 3.2
  ) +
  scale_fill_manual(values = component_colors, name = "Variance Component") +
  labs(
    x = NULL,
    y = "Adjusted R² (Variance explained)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(face = "bold"),
    strip.text.y.left = element_text(angle = 90, face = "bold"),
    plot.margin = ggplot2::margin(0.2, 0.2, 0.2, 0.2, unit = "lines")
  )



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------
# Moving window varpart - Continuous variance decomposition
# ----------------------------------------------------------

# --- Predictor variables ---
selected_variables <- c("NAO_Median_Value", "estimate")

# --- Moving window parameters ---
window_size <- 300   # years
step_size <- 30      # sliding step in years

# --- Define window sequence ---
min_year <- floor(min(joined_df_30yr$age_ce, na.rm = TRUE))
max_year <- ceiling(max(joined_df_30yr$age_ce, na.rm = TRUE))
window_starts <- seq(min_year, max_year - window_size, by = step_size)

# --- Function to compute Pure and Shared variance per window ---
calculate_pure_shared <- function(data, window_start) {
  window_end <- window_start + window_size
  subset <- data %>% filter(age_ce >= window_start, age_ce < window_end)
  
  comm <- subset %>% select(high_profile:predator) %>% filter(complete.cases(.))
  env  <- subset %>% select(all_of(selected_variables)) %>% filter(complete.cases(.))
  n <- min(nrow(comm), nrow(env))
  
  # Skip windows with insufficient data
  if (n < 5) return(NULL)
  
  comm <- comm[1:n, ]
  env  <- env[1:n, ]
  
  # Compute pure effects
  pure_effects <- map_dfr(selected_variables, function(var) {
    others <- setdiff(selected_variables, var)
    formula <- as.formula(
      paste0("comm ~ ", var, " + Condition(", paste(others, collapse = " + "), ")")
    )
    mod <- rda(formula = formula, data = env)
    tibble(
      window_start = window_start,
      window_end = window_end,
      component = ifelse(var == "estimate", "Pure Vegetation", "Pure Climate"),
      value = max(0, RsquareAdj(mod)$adj.r.squared)
    )
  })
  
  # Total explained variance
  mod_full <- rda(comm ~ ., data = env)
  total_r2 <- max(0, RsquareAdj(mod_full)$adj.r.squared)
  shared_val <- total_r2 - sum(pure_effects$value)
  
  shared <- tibble(
    window_start = window_start,
    window_end = window_end,
    component = "Shared",
    value = max(0, shared_val)
  )
  
  bind_rows(pure_effects, shared)
}

# --- Run moving window analysis ---
r2_windowed <- map_dfr(window_starts, ~calculate_pure_shared(joined_df_30yr, .x))

# --- Create midpoint labels ---
r2_windowed <- r2_windowed %>%
  mutate(
    midpoint = (window_start + window_end) / 2
  )

# --- Plot moving window variance partitioning ---
varpart_moving <- ggplot(r2_windowed, aes(x = midpoint, y = value, fill = component)) +
  geom_col(position = "stack", width = step_size) +
  scale_fill_manual(values = component_colors, name = "Variance Component") +
  labs(
    x = "Year (CE)",
    y = "Adjusted R² (Variance explained)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(face = "bold"),
    plot.margin = ggplot2::margin(0.2, 0.2, 0.2, 0.2, unit = "lines")
  ) +
  scale_x_continuous(
    breaks = seq(
      from = floor(min(r2_windowed$midpoint, na.rm = TRUE) / 100) * 100,
      to   = ceiling(max(r2_windowed$midpoint, na.rm = TRUE) / 100) * 100,
      by = 200
    )
  )




## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------
# Moving window varpart - Effect sizes
# ----------------------------------------------------------

# --- Predictor variables ---
selected_variables <- c("NAO_Median_Value", "estimate")

# --- Moving window parameters ---
window_size <- 300   # years
step_size <- 30      # sliding step in years

# --- Define window sequence ---
min_year <- floor(min(joined_df_30yr$age_ce, na.rm = TRUE))
max_year <- ceiling(max(joined_df_30yr$age_ce, na.rm = TRUE))
window_starts <- seq(min_year, max_year - window_size, by = step_size)

# --- Function to compute Pure and Shared variance per window ---
calculate_pure_shared <- function(data, window_start) {
  window_end <- window_start + window_size
  subset <- data %>% filter(age_ce >= window_start, age_ce < window_end)
  
  comm <- subset %>% select(high_profile:predator) %>% filter(complete.cases(.))
  env  <- subset %>% select(all_of(selected_variables)) %>% filter(complete.cases(.))
  n <- min(nrow(comm), nrow(env))
  
  if (n < 5) return(NULL)
  
  comm <- comm[1:n, ]
  env  <- env[1:n, ]
  
  # Compute pure effects
  pure_effects <- map_dfr(selected_variables, function(var) {
    others <- setdiff(selected_variables, var)
    formula <- as.formula(
      paste0("comm ~ ", var, " + Condition(", paste(others, collapse = " + "), ")")
    )
    mod <- rda(formula = formula, data = env)
    tibble(
      window_start = window_start,
      window_end = window_end,
      component = ifelse(var == "estimate", "Pure Vegetation", "Pure Climate"),
      value = max(0, RsquareAdj(mod)$adj.r.squared)
    )
  })
  
  # Total explained variance
  mod_full <- rda(comm ~ ., data = env)
  total_r2 <- max(0, RsquareAdj(mod_full)$adj.r.squared)
  shared_val <- total_r2 - sum(pure_effects$value)
  
  shared <- tibble(
    window_start = window_start,
    window_end = window_end,
    component = "Shared",
    value = max(0, shared_val)
  )
  
  bind_rows(pure_effects, shared)
}

# --- Run moving window analysis ---
r2_windowed <- map_dfr(window_starts, ~calculate_pure_shared(joined_df_30yr, .x))

# --- Add midpoint and split Pure Vegetation effect ---
r2_windowed <- r2_windowed %>%
  mutate(
    midpoint = (window_start + window_end) / 2,
    veg_period = case_when(
      component == "Pure Vegetation" & midpoint < 750 ~ "Pre-Coloniz. VegChange",
      component == "Pure Vegetation" & midpoint >= 750 ~ "Post-Coloniz. VegChange",
      TRUE ~ component
    )
  )

# --- Custom color palette ---
component_colors_split <- c(
  "Pre-Coloniz. VegChange" = "#054A29",       # darker green
  "Post-Coloniz. VegChange" = "#97D8C4",      # original muted green
  "Pure Climate" = "#F4B942",
  "Shared" = "#4059AD"
)

# --- Calculate effect size contrast (Vegetation - Climate) ---
effect_diff <- r2_windowed %>%
  filter(component %in% c("Pure Vegetation", "Pure Climate")) %>%
  select(midpoint, component, value) %>%
  pivot_wider(names_from = component, values_from = value) %>%
  mutate(
    effect_diff = `Pure Vegetation` - `Pure Climate`
  )




## ----fig.asp=0.9, fig.width=12------------------------------------------------------------------------------------------------------------------------------------------------------------------

# --------------------------------------------
# Define custom plot titles with island abbreviations
# --------------------------------------------
custom_titles <- c(
  "Azul" = "Azul [SM]",
  "Caldeirao" = "Caldeirão [CR]",
  "Caveiro" = "Caveiro [PI]",
  "Empadadas Norte" = "Empada. Norte [SM]",
  "Funda" = "Funda [FL]",
  "Ginjal" = "Ginjal [TE]",
  "Peixinho" = "Peixinho [PI]",
  "Prata" = "Prata [SM]",
  "Santiago" = "Santiago [SM]"
)

# --------------------------------------------
# Replace plot titles in list with custom labels
# --------------------------------------------
ls_plot_df_dca_multi_all_named <- imap(ls_plot_df_dca_multi_all, function(plot, name) {
  if (name %in% names(custom_titles)) {
    plot +
      ggtitle(custom_titles[[name]]) +
      theme(
        plot.title = element_text(size = 8, face = "bold"),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        text = element_text(size = 8)  # general fallback
      )
  } else {
    plot
  }
})


# --------------------------------------------
# Remove Furnas and Fogo plots and convert remaining plots to grobs
# --------------------------------------------
scatterplots <- 
  ls_plot_df_dca_multi_all_named[!names(ls_plot_df_dca_multi_all_named) %in% c("Furnas", "Fogo")] %>% 
  map(~ggplotGrob(.x + theme(legend.position = "none")))

# --------------------------------------------
# Convert lake coordinates to spatial points
# --------------------------------------------
lake_centers <- data.frame(island = lake_metadata$lake, lat = lake_metadata$lat, lon = lake_metadata$long)
island_points <- st_as_sf(lake_centers, coords = c("lon", "lat"), crs = 4326)

# --------------------------------------------
# Define bounding box coordinates for map panels
# --------------------------------------------
xmin <- -32.0 - 0.75 * abs(-32.0 + 24.5)
xmax <- -24.5 + 0.75 * abs(-32.0 + 24.5)
ymin <- 36.5 - 0.5 * abs(40.0 - 36.5)
ymax <- 40.0 + 0.5 * abs(40.0 - 36.5)

# Add optional buffers (set to zero here)
buffer_x <- 0 
buffer_y <- 0
xmin <- xmin - buffer_x; xmax <- xmax + buffer_x
ymin <- ymin - buffer_y; ymax <- ymax + buffer_y

# --------------------------------------------
# Create bounding boxes as spatial features
# --------------------------------------------
north_atlantic_bbox <- st_as_sfc(st_bbox(c(xmin = -100, xmax = 10, ymin = 20, ymax = 60), crs = 4326))
azores_bbox <- st_as_sfc(st_bbox(c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), crs = 4326))

# --------------------------------------------
# Load and crop coastline shapefile
# --------------------------------------------
coastline_shapefile <- st_read("data/ne_10m_coastline/ne_10m_coastline.shp")
azores_coastline <- st_crop(coastline_shapefile, azores_bbox)
north_atlantic_coastline <- st_crop(coastline_shapefile, north_atlantic_bbox)

# --------------------------------------------
# Panel A: Contextual map showing Azores within North Atlantic
# --------------------------------------------
panel_a <- ggplot() +
  geom_sf(data = north_atlantic_coastline, fill = "lightgrey", color = "black") +
  geom_sf(data = azores_bbox, fill = NA, color = "red", linetype = "dashed", size = 0.7) +
  coord_sf(xlim = c(-80, 10), ylim = c(30, 50), datum = sf::st_crs(4326), expand = FALSE) +
  theme_void() +
  ggtitle("(a) North Atlantic context of the Azores Archipelago") +
  theme(plot.title = element_text(size = 8, face = "bold")) +
  annotation_scale(location = "bl", width_hint = 0.3)

# --------------------------------------------
# Define where to position the embedded scatterplots in Panel B
# --------------------------------------------
scatterplot_positions <- list(
  c(xmax - 1.5, ymax - 1),                # Azul 
  c(xmin + 1.5, ymax - 1),                # Caldeirão
  c(xmin + 1.5, ymin + 1),                # Caveiro
  c(xmax - 1.5, (ymax + ymin) / 2),       # Empadadas Norte
  c(xmin + 1.5, (ymax + ymin) / 2),       # Funda
  c(xmin + (xmax - xmin) / 2, ymin + 1),  # Ginjal
  c(((xmin + 1.5)+(xmin + (xmax - xmin) / 2))/2, ymin + 1),     # Peixinho
  c(((xmin + (xmax - xmin) / 2)+(xmax - 1.5))/2, ymin + 1),     # Prata
  c(xmax - 1.5, ymin + 1)                 # Santiago
)

# --------------------------------------------
# Create data frame of segment lines between scatterplots and lake points
# --------------------------------------------
segment_data <- data.frame(
  x = sapply(scatterplot_positions, function(pos) pos[1]),
  y = sapply(scatterplot_positions, function(pos) pos[2]),
  xend = st_coordinates(island_points$geometry)[, 1],
  yend = st_coordinates(island_points$geometry)[, 2]
)

# --------------------------------------------
# Panel B: Zoomed-in map of Azores with embedded DCA plots
# --------------------------------------------
bg_colors <- viridis::viridis(9, option = "plasma", alpha = 0.7)

panel_b <- ggplot() +
  geom_sf(data = azores_coastline, fill = "lightblue", color = "black") +
  geom_sf(data = island_points, aes(geometry = geometry), color = "white", 
          fill = "grey", shape = 21, size = 4) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), datum = sf::st_crs(4326), expand = FALSE,
           label_graticule = "both") +
  theme_void() +
  theme(
    panel.grid.major = element_line(color = "grey70", linetype = "dotted"),
    plot.margin = unit(c(0, 0, 1, 0), "lines")
  ) +
  annotation_north_arrow(
    location = "bl", which_north = "true", 
    pad_x = unit(0.75, "in"), pad_y = unit(0.75, "in"),
    style = north_arrow_minimal()
  ) +
  geom_segment(
    data = segment_data, 
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "grey", size = 0.5
  ) +
  ggtitle("(b) Detailed map of the Azores Archipelago") +
  theme(plot.title = element_text(size = 8, face = "bold"))

# --------------------------------------------
# Overlay scatterplots at corresponding positions
# --------------------------------------------
scatterplot_scale_x <- 2 * 0.85
scatterplot_scale_y <- 1.5 * 0.75

for (i in seq_along(scatterplots)) {
  panel_b <- panel_b + 
    annotation_custom(
      scatterplots[[i]], 
      xmin = scatterplot_positions[[i]][1] - scatterplot_scale_x, 
      xmax = scatterplot_positions[[i]][1] + scatterplot_scale_x, 
      ymin = scatterplot_positions[[i]][2] - scatterplot_scale_y, 
      ymax = scatterplot_positions[[i]][2] + scatterplot_scale_y
    )
}

# --------------------------------------------
# Adjust margins and borders for layout
# --------------------------------------------
panel_a <- panel_a +
  theme(
    plot.margin = unit(c(1, 0, 0, 0), "lines"),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

panel_b <- panel_b +
  theme(plot.margin = unit(c(0, 0, 1, 0), "lines"))

# --------------------------------------------
# Combine both panels into a single composite figure
# --------------------------------------------
combined_plot_map <- panel_a / panel_b + 
  patchwork::plot_layout(heights = c(0.4, 0.6))  # Adjust vertical split

# Display the final figure
print(combined_plot_map)

# --------------------------------------------
# Export figure to file
# --------------------------------------------
ggsave(
  "figures/FIGURE_1_azores_map_DCA_all_types_two_panels.pdf",
  plot = combined_plot_map,
  dpi = 300, width = 180, height = 190, units = "mm"
)




## ----fig.asp=0.75, fig.width=12-----------------------------------------------------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------
# Combine smoothed GAM outputs and DCA scores across groups
# ----------------------------------------------------------

smooth_all_in_one <- 
  bind_rows(
    cbind(type = rep("diatoms", nrow(smooth_diatoms)), smooth_diatoms),
    cbind(type = rep("chiro", nrow(smooth_diatoms)), smooth_chiro)
  )

df_all_in_one <- 
  bind_rows(
    cbind(type = rep("diatoms", nrow(b_ls_dca_diatoms)), b_ls_dca_diatoms),
    cbind(type = rep("chiro", nrow(b_ls_dca_chiro)), b_ls_dca_chiro)
  )

# ----------------------------------------------------------
# Plot smooth trends of DCA1 scores for all groups
# ----------------------------------------------------------

plot_time.all_in_one <- 
  smooth_all_in_one %>%
  filter(.smooth == "s(age_ce)") %>%
  ggplot(aes(colour = type), group = type) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, x = age_ce, fill = type),
              alpha = 0.2, linewidth = 0) +
  geom_line(aes(x = age_ce, y = .estimate, colour = type), linewidth = 3, alpha = 1) +
  labs(y = "DCA1", title = "Chiro s(age_ce) GI model", x = "Year (CE)") +  # Added x-axis label
  scale_colour_viridis_d(name = "Groups", labels = c("Producers", "Consumers")) +
  scale_fill_viridis_d(name = "Groups", labels = c("Producers", "Consumers")) +
  theme_classic(base_size = 12) +
  lims(x = c(0, 2020), y = c(-2, 2)) +
  scale_x_continuous(n.breaks = 5, limits = c(0, 2020)) +
  theme_bw() +
  theme(
    plot.background = element_blank(),
    plot.title = element_blank(),
    legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "right",
    axis.title = element_text(size = 15),   # Increase axis titles
    axis.text = element_text(size = 13),                   # Increase axis labels
    axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
    panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )


# --------------------------------------------
# Export figure to file
# --------------------------------------------
ggsave(
  "figures/FIGURE_1b_regional_gam_DCA_all_groups.svg",
  plot = plot_time.all_in_one,
  dpi = 600, width = 12, height = 6
)




## ----fig.width = 14, fig.asp = 1.2--------------------------------------------------------------------------------------------------------------------------------------------------------------

# --------------------------------------------
# Global theme parameters (adjust once here)
# --------------------------------------------
axis_title_size <- 15
axis_text_size <- 13
plot_title_size <- 15

# --------------------------------------------
# Prepare Panel A: Ordination plot (removes legend for cleaner layout)
# --------------------------------------------
ordinations_fclusts_1 <- ordinations_fclusts_1 + 
  guides(fill = "none") +
  ggtitle("(a) Functional Structure Space") +
  theme(
    plot.margin = unit(c(1, 2, 1, 2), "pt"),
    axis.text = element_text(size = axis_text_size),
    axis.title = element_text(size = axis_title_size),
    plot.title = element_text(size = plot_title_size, face = "bold")
  )

# --------------------------------------------
# Prepare Panel B: Classification tree from PNG image
# --------------------------------------------
tree_grob <- grid::rasterGrob(png::readPNG("outputs/decision_tree_paleo.png"))
plot_decision_tree_paleo <- ggplot() +
  annotation_custom(tree_grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_void() +
  ggtitle("(b) Classification Tree") +
  theme(plot.title = element_text(size = plot_title_size, face = "bold"))

# --------------------------------------------
# Prepare Panel C: Boxplots of trophic guilds across community structures
# --------------------------------------------
facet_colors <- c(
  "CTS1" = viridis::viridis(5, direction = -1)[1],
  "CTS2" = viridis::viridis(5, direction = -1)[2],
  "CTS3" = viridis::viridis(5, direction = -1)[3],
  "CTS4" = viridis::viridis(5, direction = -1)[4],
  "CTS5" = viridis::viridis(5, direction = -1)[5]
)

boxplots <- 
  bind_cols(global_matrix_merged_paleo, select(norm_abund, c(
    euplanktonic, high_profile, low_profile, motile, algivore, 
    detritivore, plantivore, predator))) %>%
  pivot_longer(
    cols = c(euplanktonic, high_profile, low_profile, motile, algivore, 
             detritivore, plantivore, predator),
    names_to = "name",
    values_to = "value"
  ) %>%
  mutate(
    amd_clusts = as.factor(amd_clusts),
    name = fct_relevel(name,
      "euplanktonic", "high_profile", "low_profile", "motile",
      "algivore", "detritivore", "plantivore", "predator"
    )
  ) %>%
  ggplot(aes(name, value)) +
  geom_boxplot() +
  scale_x_discrete(labels = c(
    "euplanktonic" = "EU", "high_profile" = "HP", "low_profile" = "LP", "motile" = "MOT",
    "algivore" = "AL", "detritivore" = "DET", "plantivore" = "PL", "predator" = "PR"
  )) +
  theme_bw() +
  labs(x = "Community Trophic Guilds", y = "Relative Abundance") +
  guides(fill = "none") +
  facet_wrap2(
    ncol = 6,
    ~ factor(amd_clusts, labels = c("CTS1", "CTS2", "CTS3", "CTS4", "CTS5")),
    strip = strip_themed(
      text_x = element_text(color = "black"),
      background_x = list(
        "CTS1" = element_rect(fill = scales::alpha(facet_colors["CTS1"], 0.8)),
        "CTS2" = element_rect(fill = scales::alpha(facet_colors["CTS2"], 0.8)),
        "CTS3" = element_rect(fill = scales::alpha(facet_colors["CTS3"], 0.8)),
        "CTS4" = element_rect(fill = scales::alpha(facet_colors["CTS4"], 0.8)),
        "CTS5" = element_rect(fill = scales::alpha(facet_colors["CTS5"], 0.8))
      )
    )
  ) +
  ggtitle("(c) Guild Composition by Structure") +
  theme(
    plot.margin = unit(c(1, 2, 1, 2), "pt"),
    plot.title = element_text(size = plot_title_size, face = "bold"),
    axis.text = element_text(size = axis_text_size),
    axis.title = element_text(size = axis_title_size),
    legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
    legend.position = "none",
    panel.grid = element_blank(),
    axis.ticks = element_line(colour = "grey70", linewidth = 0.2)
  )

# --------------------------------------
# Prepare Panel D: Vegetation change with GAM smooth
# --------------------------------------
vegetation_smooth_plot <- ggplot(interpolated_df, aes(x = age_ce, y = estimate)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "steelblue", alpha = 0.2) +
  theme_minimal(base_size = 14) +
  labs(title = "Vegetation Change Regional Smooth", x = "Year CE", y = "Estimate")

df_p1_smooth <- ggplot_build(vegetation_smooth_plot)$data[[1]] %>%
  filter(!is.na(y)) %>%
  mutate(y_rescaled = scales::rescale(y, to = c(0, 1)))

df_p1_ribbon <- ggplot_build(vegetation_smooth_plot)$data[[2]] %>%
  filter(!is.na(ymin) & !is.na(ymax)) %>%
  mutate(
    ymin_rescaled = scales::rescale(ymin, to = c(0, 1)),
    ymax_rescaled = scales::rescale(ymax, to = c(0, 1))
  )

orig_range <- range(df_p1_smooth$y, na.rm = TRUE)

p_combined <- plot_stacked_fclusts_time_lines +
  geom_line(data = df_p1_smooth, aes(x = x, y = y_rescaled),
            inherit.aes = FALSE, color = "white", size = 2) +
  scale_y_continuous(
    name = "Proportion of Community Trophic Structures",
    sec.axis = sec_axis(
      trans = ~ scales::rescale(., from = c(0, 1), to = orig_range),
      name = "Vegetation Change (GAM smooth - Pollen signal)"
    )
  ) +
  theme_minimal(base_size = 12) +
  ggtitle("(d) Proportion of Community Trophic Structures Over Time") +
  theme(
    plot.title = element_text(size = plot_title_size, face = "bold"),
    axis.text = element_text(size = axis_text_size),
    axis.title = element_text(size = axis_title_size)
  )

# --------------------------------------
# Prepare Panel E: Richness boxplots
# --------------------------------------
box_plot_clusts_mod <- box_plot_clusts_mod +
  scale_x_discrete(labels = c("CTS5", "CTS4", "CTS3", "CTS2", "CTS1")) +
  ggtitle("(e) Number of Species per Community Trophic Structure") +
  theme(
    plot.title = element_text(size = plot_title_size, face = "bold"),
    axis.text = element_text(size = axis_text_size),
    axis.title = element_text(size = axis_title_size)
  )

# --------------------------------------------
# Define patchwork layout
# --------------------------------------------
layout <- "
AABB
AABB
CCCC
CCCC
"

# --------------------------------------------
# Assemble Panels A–C
# --------------------------------------------
figure_1_ord_rf_time_clusts <- (
  ordinations_fclusts_1 +
  plot_decision_tree_paleo +
  boxplots
) +
  patchwork::plot_layout(design = layout, heights = c(40, 40, 30), guides = "collect") +
  patchwork::plot_annotation() &
  theme(legend.position = 'none')

# --------------------------------------------
# Final assembly of all five panels
# --------------------------------------------
final_figure <- (
  figure_1_ord_rf_time_clusts /
  (
    (p_combined + theme(legend.position = "none")) +
    (box_plot_clusts_mod + theme(legend.position = "none"))
  ) +
  plot_layout(guides = "collect", heights = c(8, 5)) +
  plot_annotation(tag_levels = NULL)
)
print(final_figure)
# --------------------------------------------
# Save to file
# --------------------------------------------
ggsave(plot = final_figure, 
       filename = "figures/FIGURE_2_composite_figure.png", 
       dpi = 600,
       width = 14)




## ----fig.asp = 0.5, fig.width = 10--------------------------------------------------------------------------------------------------------------------------------------------------------------

# --------------------------------------------
# Set common visual parameters (edit these to update all plots)
# --------------------------------------------
axis_title_size <- 13
axis_text_size  <- 11
legend_text_size <- 10
plot_title_size <- 14

# --------------------------------------------
# Helper: Standardize and summarize species richness
# --------------------------------------------
summarize_richness <- function(df, producer = TRUE) {
  df %>%
    select(-name) %>%
    group_by(lake, core_depth_id, fgroup) %>%
    summarise(
      species_richness = sum(value > 0, na.rm = TRUE),
      total_abundance  = sum(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(lake, fgroup) %>%
    mutate(
      standardized_species_richness = (species_richness - min(species_richness, na.rm = TRUE)) / 
                                      (max(species_richness, na.rm = TRUE) - min(species_richness, na.rm = TRUE))
    ) %>%
    left_join(df_age_model, by = "core_depth_id") %>%
    filter(if (producer) !fgroup %in% c("algivore", "detritivore", "plantivore", "predator") else fgroup %in% c("algivore", "detritivore", "plantivore", "predator"))
}

# --------------------------------------------
# Global consumer richness trend
# --------------------------------------------
cons_glob_div_plt <- summarize_richness(df_comm, producer = FALSE) %>%
  ggplot(aes(age_ce, standardized_species_richness)) +
  geom_smooth(aes(color = fgroup), se = TRUE) +
  geom_smooth(aes(linetype = "Overall mean"), color = "red", alpha = 0.7, se = FALSE) +
  scale_color_viridis_d(name = "Consumer guilds", option = "B", direction = -1, alpha = 0.5) +
  scale_linetype_manual(name = "Overall mean", values = c("Overall mean" = 3)) +
  labs(x = "Age (CE)", y = "Standardized Species Richness (0–1)") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = axis_title_size),
    axis.text  = element_text(size = axis_text_size),
    plot.title = element_text(size = plot_title_size, face = "bold")
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  lims(x = c(0, 2010))

# --------------------------------------------
# Per-lake consumer richness trends
# --------------------------------------------
cons_lake_div_plt <- summarize_richness(df_comm, producer = FALSE) %>%
  ggplot(aes(age_ce, standardized_species_richness)) +
  facet_wrap(~lake, scales = "free", nrow = 3) +
  geom_point(aes(color = fgroup), alpha = 0.6) +
  geom_smooth(aes(color = fgroup), se = FALSE) +
  geom_smooth(aes(linetype = "Overall mean"), color = "red", alpha = 0.7, se = FALSE) +
  scale_color_viridis_d(name = "Consumer guilds", option = "B", direction = -1, alpha = 0.5) +
  scale_linetype_manual(name = "Overall mean", values = c("Overall mean" = 3)) +
  labs(x = "Age (CE)", y = "Standardized Species Richness (0–1)") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = axis_title_size),
    axis.text  = element_text(size = axis_text_size),
    plot.title = element_text(size = plot_title_size, face = "bold")
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  lims(x = c(0, 2010))

# --------------------------------------------
# Global producer richness trend
# --------------------------------------------
prod_glob_div_plt <- summarize_richness(df_comm, producer = TRUE) %>%
  ggplot(aes(age_ce, standardized_species_richness)) +
  geom_smooth(aes(color = fgroup), se = TRUE) +
  geom_smooth(aes(linetype = "Overall mean"), color = "red", alpha = 0.7, se = FALSE) +
  scale_color_viridis_d(name = "Producer guilds", option = "D", direction = -1, alpha = 0.5) +
  scale_linetype_manual(name = "Overall mean", values = c("Overall mean" = 3)) +
  labs(x = "Age (CE)", y = "Standardized Species Richness (0–1)") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = axis_title_size),
    axis.text  = element_text(size = axis_text_size),
    plot.title = element_text(size = plot_title_size, face = "bold")
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  lims(x = c(0, 2010))

# --------------------------------------------
# Per-lake producer richness trends
# --------------------------------------------
prod_lake_div_plt <- summarize_richness(df_comm, producer = TRUE) %>%
  ggplot(aes(age_ce, standardized_species_richness)) +
  facet_wrap(~lake, scales = "free", nrow = 3) +
  geom_point(aes(color = fgroup), alpha = 0.6) +
  geom_smooth(aes(color = fgroup), se = FALSE) +
  geom_smooth(aes(linetype = "Overall mean"), color = "red", alpha = 0.7, se = FALSE) +
  scale_color_viridis_d(name = "Producer guilds", option = "D", direction = -1, alpha = 0.5) +
  scale_linetype_manual(name = "Overall mean", values = c("Overall mean" = 3)) +
  labs(x = "Age (CE)", y = "Standardized Species Richness (0–1)") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = axis_title_size),
    axis.text  = element_text(size = axis_text_size),
    plot.title = element_text(size = plot_title_size, face = "bold")
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  lims(x = c(0, 2010))

# --------------------------------------------
# Optional: Remove y-axis for consumer global plot to align with producer
# --------------------------------------------
cons_glob_div_plt <- cons_glob_div_plt +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank()
  )

# --------------------------------------------
# Combine producer and consumer global plots
# --------------------------------------------
plot_div_per_fgroup_regional_scale <- 
  prod_glob_div_plt + cons_glob_div_plt +
  patchwork::plot_layout(guides = "collect", widths = c(1, 1)) +
  patchwork::plot_annotation(
    tag_levels = 'A',
    theme = theme(
      legend.position = "bottom",
      legend.text = element_text(size = legend_text_size)
    )
  )

# --------------------------------------------
# Save output
# --------------------------------------------
ggsave(
  plot = plot_div_per_fgroup_regional_scale,
  filename = "figures/FIGURE_3_plot_div_per_fgroup_regional_scale.png",
  dpi = 600,
  width = 10, height = 6
)



## ----fig.width = 10, fig.asp = 1----------------------------------------------------------------------------------------------------------------------------------------------------------------
# STEP 1: Define historical transition years
phase_lines <- c(750, 1050, 1450, 1750)

# STEP 2: Plot A – Moving window R² (unchanged from earlier setup)
varpart_moving <- ggplot(r2_windowed, aes(x = midpoint, y = value, fill = veg_period)) +
  geom_col(position = "stack", width = step_size) +
  scale_fill_manual(values = component_colors_split, name = "Variance Component") +
  scale_x_continuous(
    breaks = seq(
      from = floor(min(r2_windowed$midpoint) / 100) * 100,
      to   = ceiling(max(r2_windowed$midpoint) / 100) * 100,
      by = 200
    )
  ) +
  geom_vline(xintercept = phase_lines, linetype = "dotted", color = "gray30") +
  labs(
    x = NULL,
    y = "Adjusted R² (Variance explained)",
    title = "(a) Variance partitioning over time"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = ggplot2::margin(0.2, 0.2, 0, 0.2, unit = "lines")
  )

# STEP 3: Compute effect size difference + color classification
effect_diff <- r2_windowed %>%
  filter(component %in% c("Pure Vegetation", "Pure Climate")) %>%
  select(midpoint, component, value) %>%
  pivot_wider(names_from = component, values_from = value) %>%
  mutate(
    effect_diff = `Pure Vegetation` - `Pure Climate`,
    effect_group = case_when(
      effect_diff > 0 & midpoint < 750 ~ "VegChange > Climate (pre-750)",
      effect_diff > 0 & midpoint >= 750 ~ "VegChange > Climate (post-750)",
      effect_diff < 0 ~ "Climate > VegChange",
      TRUE ~ "Equal"
    )
  )

# STEP 4: Define custom color palette
effect_colors_split <- c(
  "VegChange > Climate (pre-750)" = "#054A29",    # dark green
  "VegChange > Climate (post-750)" = "#6DC6B6",   # lighter green
  "Climate > VegChange" = "#F4B942"               # orang
)

# STEP 5: Plot B – Effect size plot with split greens
effect_diff_plot <- ggplot(effect_diff, aes(x = midpoint, y = effect_diff, fill = effect_group)) +
  geom_col(width = step_size) +
  scale_fill_manual(values = effect_colors_split, name = "Dominant Driver") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = phase_lines, linetype = "dotted", color = "gray30") +
  scale_x_continuous(
    breaks = seq(
      from = floor(min(effect_diff$midpoint) / 100) * 100,
      to   = ceiling(max(effect_diff$midpoint) / 100) * 100,
      by = 200
    )
  ) +
  labs(
    x = "Year (CE)",
    y = "Effect size (Vegetation − Climate)",
    title = "(b) Relative dominance of VegChange vs. climate"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = ggplot2::margin(0.2, 0.2, 0.2, 0.2, unit = "lines")
  )


# Apply title formatting in each subplot
original_varpart <- original_varpart + 
  ggtitle("(a) Historical phases") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0))

varpart_moving <- varpart_moving +
  ggtitle("(b) Moving window") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0))

effect_diff_plot <- effect_diff_plot +
  ggtitle("(c) Climate vs VegChange") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0))

# Build the combined plot
combined_varplot <- (
  (original_varpart + ggtitle("(a) Historical phases")) /
  (varpart_moving + ggtitle("(b) Moving window")) /
  (effect_diff_plot + ggtitle("(c) Climate vs VegChange"))
) +
  #plot_layout(guides = "collect") &
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0),
    plot.margin = ggplot2::margin(6, 6, 6, 6)
  )

# Print the plot
print(combined_varplot)

# --------------------------------------------
# Save output
# --------------------------------------------
ggsave(
  plot = combined_varplot,
  filename = "figures/FIGURE_4_multivariate_analysis.png",
  dpi = 600,
  width = 14, height = 10
)



## ----fig.width=12, fig.asp=0.5------------------------------------------------------------------------------------------------------------------------------------------------------------------

# --------------------------------------
# Plot 1: NH Summer Temperature (Büntgen tree-ring reconstruction)
# --------------------------------------
nhst_buntgen_plot <- ggplot(nhst_buntgen, aes(x = age_ce, y = Rmean)) +
  geom_line(size = 1, color = "grey") +  # Main raw line
  geom_point(aes(color = sign), size = 1) +  # Positive (red) / negative (blue) deviations
  scale_color_manual(values = c("positive" = "red", "negative" = "blue")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +  # zero reference

  # Add 50-year smoothed trend line
  geom_line(data = nhst_buntgen_binned, aes(x = bin_50yr, y = avg_value),
            color = "steelblue", linewidth = 2) +

  theme_minimal(base_size = 14) +
  labs(
    title = "NH Summer Temperature reconstruction (Tree Rings)",
    x = "Year CE",
    y = "Temperature Anomaly",
    color = NULL
  ) +
  theme(legend.position = "none")

# --------------------------------------
# Plot 2: North Atlantic Oscillation (NAO) median values
# --------------------------------------
nao_plot <- ggplot(nao_hernandez, aes(x = age_ce, y = NAO_Median_Value)) +
  geom_line(size = 1, color = "grey") +  # Raw NAO data
  geom_point(aes(color = sign), size = 1) +  # Red/blue anomalies
  scale_color_manual(values = c("positive" = "red", "negative" = "blue")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +

  # Add 50-year smoothed trend line
  geom_line(data = nao_binned, aes(x = bin_50yr, y = avg_value),
            color = "steelblue", linewidth = 2) +

  theme_minimal(base_size = 14) +
  labs(
    title = "NAO Anomaly",
    x = "Year CE",
    y = "NAO Median",
    color = NULL
  ) +
  theme(legend.position = "none")

# --------------------------------------
# Custom Theme for Plot Tag Styling
# --------------------------------------
tag_theme <- theme(
  plot.tag.position = c(0, 1),
  plot.tag = element_text(
    hjust = 0,
    vjust = -2,
    face = "bold",
    size = 11,
    margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 10, unit = "pt")  # Left indent for tag
  ),
  plot.margin = ggplot2::margin(20, 10, 10, 10, unit = "pt")  # Standard margins
)

# --------------------------------------
# Assemble Tagged and Themed Panels
# --------------------------------------
nhst_buntgen_panel <- nhst_buntgen_plot + 
  labs(tag = "(a) NH Summer Temperature Anomaly reconstruction (Tree Rings)", title = NULL) + 
  tag_theme +
  theme(axis.title.x = element_blank())  # remove x-axis label

nao_panel <- nao_plot + 
  labs(tag = "(b) NAO Index", title = NULL) + 
  tag_theme +
  theme(axis.title.x = element_blank())  # remove x-axis label

# --------------------------------------
# Combine Panels Vertically with patchwork
# --------------------------------------
final_plot_2 <- 
  (nhst_buntgen_panel + nao_panel) +
  patchwork::plot_layout(ncol = 1, nrow = 3,
                         heights = c(20, 20, 20)) +  # Equal height for each plot
  patchwork::plot_annotation(
    tag_levels = NULL,
    theme = theme(
      plot.margin = ggplot2::margin(10, 10, 30, 10, unit = "pt"),  # extra bottom space
      plot.caption = element_text(hjust = 0.5, face = "bold", size = 12)
    )
  )

# --------------------------------------
# Export to PNG
# --------------------------------------

# Set PNG output for saving the tree plot
png("supplementary/supplementary_figure_2.png", width = 8, height = 8, units = "in", res = 600)

# Plot
final_plot_2

# Close the graphics device
dev.off()




## ----fig.width = 8, fig.asp = 1.5---------------------------------------------------------------------------------------------------------------------------------------------------------------

# --------------------------------------
# Plot 1: Group-Specific Pollen Trends (Tree vs Herb)
# --------------------------------------
group_trend_plot <- ggplot(combined_df_pollen, aes(x = age_ce, y = percent, color = group)) +
  geom_point(alpha = 0.2, size = 1) +  # Raw points with transparency
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, k = 30),  # GAM smoothing with 30 knots
    se = TRUE,
    linewidth = 1.5
  ) +
  scale_color_manual(
    values = c("Tree" = "#3a6520", "Herb" = "#dbbc3f"),  # Custom color for each pollen group
    labels = c("Tree" = "Tree", "Herb" = "Herb"),
    name = "Pollen Group"
  ) +
  labs(
    title = "Vegetation change",
    x = "Age (CE)",
    y = "% Abundance"
  ) +
  coord_cartesian(xlim = c(0, 2010)) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )

# --------------------------------------
# Plot 2: Regional Pollen Turnover (GI model smooth)
# --------------------------------------
regional_smooth_plot <- interpolated_df %>%
  ggplot() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = age_ce),
              alpha = 0.2) +  # Confidence interval
  geom_line(aes(x = age_ce, y = estimate), lwd = 1) +  # Smooth estimate line
  labs(y = "DCA1", title = "s(age_ce) GI model") +
  scale_color_manual("Time-Interval", values = lakes_col) +  # Optional (no visible use)
  theme_minimal(base_size = 12) +
  lims(y = c(-1, 1)) +
  labs(title = NULL) +
  theme(
    panel.grid.minor.y = element_blank(),
    legend.position = "none"
  )

# --------------------------------------
# Panel A: Regional pollen turnover with label
# --------------------------------------
regional_panel_A <- regional_smooth_plot +
  labs(title = "(a) Regional Pollen Community Turnover") +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10)
  )

# --------------------------------------
# Panel B: Tree vs Herb pollen dynamics with label
# --------------------------------------
group_panel_B <- group_trend_plot +
  labs(title = "(b) Tree and Herb pollen dynamics") +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10)
  )

# --------------------------------------
# Combine Panels into Final Supplementary Figure
# --------------------------------------
final_combined_plot_pollen <- 
  regional_panel_A + 
  group_panel_B +
  patchwork::plot_layout(nrow = 2, heights = c(1, 1)) +
  patchwork::plot_annotation(
    theme = theme(
      plot.margin = ggplot2::margin(10, 10, 10, 10),
      plot.title = element_text(face = "bold", size = 14)
    )
  )

# --------------------------------------
# Export to PNG
# --------------------------------------
png("supplementary/supplementary_figure_3_vegetation_change.png", 
    width = 7.29, height = 4.51, units = "in", res = 600)

# Render the plot to file
final_combined_plot_pollen

# Close the graphics device
dev.off()



## ----fig.asp = 1.5, fig.width = 12--------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------
# Convert PNG tree plot into a ggplot-compatible object
# -----------------------------------------------
tree_grob_paleo <- grid::rasterGrob(png::readPNG("outputs/decision_tree_paleo.png"))

# Create a ggplot with the decision tree image embedded as a background
tree_image_panel <- ggplot() +
  annotation_custom(tree_grob_paleo, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_void() +
  ggtitle("Functional group classification by ecosystem type")  # Title for Panel D

# -----------------------------------------------
# Clean up lake-scale density plot (remove redundant legend)
# -----------------------------------------------
lake_scale_density_plot <- plot_density_fclusts_time_lines_1B + 
  guides(fill = "none")  # Remove fill legend (i.e., for functional groups)

# -----------------------------------------------
# Define layout structure for patchwork combination
# -----------------------------------------------
patchwork_layout_design <- "
AAAA
DDDD
DDDD
DDDD
"

# -----------------------------------------------
# Assemble Supplementary Figure 4
# -----------------------------------------------

# Panel A: Regional scale
panel_regional_stacked <- plot_stacked_fclusts_time_lines +
  ggtitle(" (a) Community Trophic Structures over time [regional scale]") +
  theme(
    plot.margin = unit(c(1, 2, 1, 2), "pt"),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold")
  )

# Panel B: Lake scale
panel_lake_density <- lake_scale_density_plot +
  ggtitle(" (b) Community Trophic Structures over time [lake scale]") +
  theme(
    plot.margin = unit(c(1, 2, 1, 2), "pt"),
    legend.position = "top",
    legend.key.size = unit(0.5, "lines"),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold")
  )
# Combine panels using patchwork with layout control
supplementary_figure_4 <- 
  (panel_regional_stacked + panel_lake_density) +
  patchwork::plot_layout(
    design = patchwork_layout_design,
    heights = c(40, 80),
    guides = "collect"
  ) +
  #patchwork::plot_annotation(tag_levels = 'A') &
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "pt"),  # Minimize outer spacing
    legend.position = 'bottom'               # Merge legends at the bottom
  )

# -----------------------------------------------
# Display the combined supplementary figure
# -----------------------------------------------
print(supplementary_figure_4)

# -----------------------------------------------
# Save to file
# -----------------------------------------------
ggsave(
  plot = supplementary_figure_4,
  filename = "supplementary/supplementary_figure_4.png",
  dpi = 600,
  width = 12,
  height = 16
)



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Set PNG output for saving the tree plot
png("supplementary/supplementary_figure_5_decision_tree.png", width = 7.29, height = 4.51, units = "in", res = 600)

# Plot variable importance from random forest
varImpPlot(RFfit_paleo, main = "")

# Close the graphics device
dev.off()



## ----fig.width = 12, fig.asp = 1----------------------------------------------------------------------------------------------------------------------------------------------------------------

# Set PNG output for saving the tree plot
png("supplementary/supplementary_figure_6_global_rda.png", width = 8, height = 8, units = "in", res = 600)

# Plot Global Rda with selected variables
global_rda

# Close the graphics device
dev.off()


