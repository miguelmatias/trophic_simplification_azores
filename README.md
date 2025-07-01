# Paleo-Trophic Structures Analysis

This repository contains scripts and data for analyzing trophic structure changes in Azorean lake ecosystems over the past two millennia. It combines paleoecological proxies, statistical modeling, and data visualization to assess changes in biodiversity and community organization.

![Community Structure](figures/FIGURE_2_composite_figure.png)

## Project Structure

```
├── data/                     # Local data files (climate reconstructions, pollen data, etc.)
├── figures/                  # Output plots and figures
├── supplementary/            # Output plots and figures
├── outputs/                  # Output intermediate figures 
├── functions/                # Custom R functions
└── main_script.Rmd           # Primary analysis script (R Markdown)
```

## Requirements

All scripts are written in **R (≥4.2.0)** and rely on the following major packages:

### 🧰 General Data Handling
- `tidyverse` (includes `dplyr`, `ggplot2`, `tidyr`, `tibble`, `purrr`, `forcats`)
- `readxl`, `WriteXLS`, `googlesheets4` — for reading Excel and Google Sheets
- `zoo`, `tidyquant` — for time series and interpolation
- `rioja` — for paleoecological utilities
- `parallel`, `furrr` — for parallel processing
- `scales`, `png`, `Hmisc`

### 📊 Visualization & Plotting
- `ggplot2` (via `tidyverse`)
- `ggrepel`, `ggpmisc`, `viridis`, `cowplot`, `patchwork`
- `ggtext`, `ggridges`, `ggh4x`, `ggspatial`, `ggordiplots`, `ggnewscale`, `gridExtra`, `grid`

### 📈 Statistical Modeling
- `mgcv`, `scam`, `gratia` — for GAMs and smoothers
- `vegan`, `codyn` — for community ecology and dynamics
- `broom` — for tidying model outputs
- `rshift` — for changepoint detection

### 🔗 Network Analysis
- `igraph`, `ggnetwork`, `qgraph`

### 🤖 Machine Learning & Tree Models
- `caret`, `randomForest`, `xgboost`, `earth`, `e1071`
- `rpart`, `rpart.plot`, `evtree` — for decision trees
- `DALEX`, `DALEXtra` — for model explainability
- `irr` — inter-rater reliability

### 🔬 Phylogenetics & Spatial
- `ape` — for phylogenetic trees
- `sf` — for spatial vector data

### 📂 Project-Specific
- Custom functions sourced from `functions/source_custom_functions.R`

### Install all packages at once

```r
install.packages(c(
  "tidyverse", "readxl", "WriteXLS", "googlesheets4", "zoo", "tidyquant", "rioja", 
  "parallel", "furrr", "scales", "png", "Hmisc", 
  "ggrepel", "ggpmisc", "viridis", "cowplot", "patchwork", "ggtext", 
  "ggridges", "ggh4x", "ggspatial", "ggordiplots", "ggnewscale", "gridExtra", "grid", 
  "vegan", "codyn", "broom", "mgcv", "scam", "gratia", "rshift",
  "igraph", "ggnetwork", "qgraph",
  "caret", "randomForest", "xgboost", "earth", "e1071", "rpart", "rpart.plot", "evtree", 
  "DALEX", "DALEXtra", "irr", 
  "ape", "sf"
))
```

---

## Data

Ensure the following datasets are available in the `data/` folder or downloaded from the corresponding sources:

- Lake metadata (`table_lake_metadata.xlsx`)
- Abundance and diversity of guilds
- NAO and SST reconstructions
- HGAM output from vegetation turnover
- Pollen counts from regional cores

## Running the analysis

Open `main_script.Rmd` in RStudio and execute the script chunk-by-chunk or use **Knit** to render a full report.

## Outputs

- Figures exported to the `figures/` folder
- Supplementary Figures exported to the `supplementary/` folder

## Author

Miguel Matias, Museo Nacional de Ciencias Naturales (MNCN-CSIC)
