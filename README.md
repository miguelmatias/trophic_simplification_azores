# Paleo-Trophic Structures Analysis

This repository contains scripts and data for analyzing trophic structure changes in Azorean lake ecosystems over the past two millennia. It combines paleoecological proxies, statistical modeling, and data visualization to assess changes in biodiversity and community organization.

![Community Structure](figures/FIGURE_2_composite_figure.png)

## Project Structure

```
├── data/                     # Local data files (climate reconstructions, pollen data, etc.)
├── figures/                  # Output plots and figures
├── supplementary/            # Output supplementary plots and figures
├── outputs/                  # Output intermediate figures 
├── functions/                # Custom R functions
└── main_script.Rmd           # Primary analysis script (R Markdown)
```

## Requirements

All scripts are written in **R (≥4.2.0)** and rely on the following major packages:

### 🧰 General Data Handling
- `tidyverse` (includes `dplyr`, `ggplot2`, `tidyr`, `tibble`, `purrr`, `forcats`)
- `zoo` — for time series manipulation and interpolation
- `parallel`, `furrr` — for parallel computing
- `scales`, `png`, `Hmisc`, `ggtext`, `scico` — for plotting utilities, formatting, and color palettes
- `conflicted` — to manage function name conflicts explicitly
- `sf` — for spatial vector data

### 📊 Visualization & Plotting
- `ggplot2` (via `tidyverse`)
- `ggrepel`, `ggpmisc`, `viridis`, `cowplot`, `patchwork`
- `ggtext`, `ggridges`, `ggh4x`, `ggspatial`, `ggordiplots`, `ggnewscale`, `gridExtra`

### 📈 Statistical Modeling
- `mgcv`, `scam`, `gratia` — for Generalized Additive Models and shape-constrained smoothers
- `vegan`, `codyn` — for ordination, diversity, and community dynamics
- `broom` — to convert model outputs into tidy data frames

### 🤖 Machine Learning & Tree Models
- `caret`, `e1071`, `evtree`, `rpart`, `randomForest`, `rpart.plot`, `earth`, `xgboost`
- `DALEX`, `DALEXtra` — for model interpretability
- `irr` — for inter-rater agreement metrics
- `ape` — for tree analysis

---

### 🔧 Setup

You can install all required packages using the following R code:

```r
required_packages <- c(
  "tidyverse", "zoo", "parallel", "furrr", "scales", "png", "ggtext", "Hmisc", "scico", "conflicted",
  "ggrepel", "ggpmisc", "viridis", "cowplot", "patchwork", "gridExtra", "ggnewscale", "ggridges",
  "ggh4x", "ggspatial", "ggordiplots", "vegan", "codyn", "broom", "mgcv", "scam", "gratia",
  "caret", "e1071", "evtree", "rpart", "randomForest", "rpart.plot", "earth", "xgboost",
  "DALEX", "DALEXtra", "irr", "ape", "sf"
)

installed <- required_packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(required_packages[!installed])
}

```
---

# Data

Ensure the following datasets are available in the `data/` directory before running the scripts:

### 🗺 Lake Metadata
- `table_lake_metadata.csv` — lake names, coordinates, area, and other physical descriptors

### 🌿 Guild Abundances & Diversity
- `table_abundance_guilds.csv` — relative abundance of trophic guilds per sample
- `table_diversity_guilds.csv` — diversity metrics (e.g., richness, evenness) per guild

### 🌡 Climate & Environmental Reconstructions
- `table_nao_hernandez.csv` — NAO index (Hernandez et al.)
- `table_nhst_buntgen.csv` — Northern Hemisphere summer temperatures (Büntgen et al.)
- `table_sst_jiang.csv` — Sea Surface Temperature reconstructions (Jiang et al.)
- `table_sst_abrantes.csv` — Sea Surface Temperature reconstructions (Abrantes et al.)

### 🍃 Vegetation & Pollen Data
- `table_pollen_azores.csv` — regional pollen group counts and composition

All data files must be placed in the `data/` subfolder of the project root directory.

## Running the analysis

Open `main_script.Rmd` in RStudio and execute the script chunk-by-chunk or use **Knit** to render a full report.

## Outputs

- Figures exported to the `figures/` folder
- Supplementary Figures exported to the `supplementary/` folder

## Author

Miguel Matias, Museo Nacional de Ciencias Naturales (MNCN-CSIC)
