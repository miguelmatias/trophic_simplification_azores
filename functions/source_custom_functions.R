
library(e1071)
library(evtree)
library(irr)
library(ggplot2)
library(rpart)
library(randomForest)
library(dplyr)
library(gridExtra)
library(grid)
library(ggridges)
library(earth)

#' Check Collinearity of Predictor Variables Using VIF
#'
#' Computes Variance Inflation Factors (VIF) for a set of explanatory variables
#' using redundancy analysis (RDA) and prints a warning if any VIF exceeds a specified threshold.
#'
#' @param vars A data frame or tibble of numeric explanatory variables to check for collinearity.
#' @param threshold Numeric. The VIF threshold above which variables are considered collinear. Default is 10.
#'
#' @return A named numeric vector of VIF values for each variable.
#' If any VIF exceeds the threshold, a warning is printed.
#'
#' @examples
#' \dontrun{
#'   df <- data.frame(a = rnorm(100), b = rnorm(100), c = rnorm(100))
#'   check_collinearity(df, threshold = 5)
#' }
#'
#' @importFrom vegan rda vif.cca
#' @export
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

#' Flip DCA1 Axis for Temporal Interpretability
#'
#' This function ensures that the DCA1 axis is temporally interpretable
#' by checking whether more recent samples (younger ages) have lower DCA1
#' values than older ones. If not, the DCA1 axis is flipped.
#'
#' @param data A data frame containing at least the columns `age_ce` and `DCA1`.
#'        `age_ce` should represent calendar age in CE (Common Era), and
#'        `DCA1` the first Detrended Correspondence Analysis axis.
#'
#' @return A data frame with DCA1 possibly flipped (multiplied by -1)
#'         to ensure interpretability in time series.
#'
#' @details The function splits the dataset into two halves based on age:
#' recent (younger than midpoint) and past (older than midpoint). If the mean
#' DCA1 in the past is greater than in the recent samples, the DCA1 axis is
#' flipped. Otherwise, it is returned unchanged.
#'
#' @export
convert_DCA <- function(data) {
  cut_off <- max(data$age_ce, na.rm = TRUE) - 
    ((max(data$age_ce, na.rm = TRUE) - min(data$age_ce, na.rm = TRUE)) / 2)
  
  # Mean DCA1 for recent and past halves of the core
  m_recent <- data[data$age_ce > cut_off & !is.na(data$age_ce), "DCA1"] %>%
    unlist() %>% mean()
  m_past <- data[data$age_ce < cut_off & !is.na(data$age_ce), "DCA1"] %>%
    unlist() %>% mean()
  
  # Flip sign of DCA1 if older samples have higher values than recent ones
  if (m_past > m_recent) {
    data <- data %>% mutate(DCA1 = DCA1 * -1)
  }
  return(data)
}


#' Subsample Time Series to Balance Timepoints Across Bins
#'
#' This function subsamples a time series to balance the number of timepoints across bins of a specified width. It randomly samples a defined number of data points from each bin.
#'
#' @param data A data frame containing the time series data.
#' @param x A character string representing the column name for the time variable (e.g., years).
#' @param y A character string representing the column name for the dependent variable (e.g., DCA axis value).
#' @param binwidth A numeric value specifying the width of the bins (e.g., 10 for 10-year bins).
#' @param n_samples_per_bin A numeric value specifying the number of samples to draw from each bin. If a bin contains fewer points than `n_samples_per_bin`, all points from that bin are returned.
#'
#' @return A data frame containing the subsampled data with balanced timepoints across bins.
#' 
#' @examples
#' # Sample data
#' data <- data.frame(
#'   x = seq(1950, 2020, length.out = 100),
#'   y = rnorm(100)
#' )
#' 
#' # Subsample with 10-year bins and 5 samples per bin
#' subsample_time_series(data = data, x = "x", y = "y", binwidth = 10, n_samples_per_bin = 5)
#'
#' @import dplyr
#' @import ggplot2
#' @export
subsample_time_series <- function(data, x, y, binwidth, n_samples_per_bin, seed = 123) {
  if (!is.null(seed)) set.seed(seed)
  
  # Create a bin column based on the binwidth
  data$bin <- cut(data[[x]], breaks = seq(min(data[[x]]), max(data[[x]]), by = binwidth), include.lowest = TRUE)
  
  # Randomly subsample the same number of points from each bin
  balanced_data <- data %>%
    group_by(bin) %>%
    sample_n(size = min(n(), n_samples_per_bin), replace = FALSE)
  
  # Return the subsampled data
  return(balanced_data)
}

#' Fix Slumps by Adding Decimals to Duplicate Ages
#'
#' Resolves duplicate \code{age_ce} values (e.g. sediment slumps) by adding
#' small decimal offsets so each row has a unique age. Used for age-depth
#' or time-series consistency.
#'
#' @param .ages A data frame with column \code{age_ce} (numeric).
#' @return The same data frame with \code{age_ce} modified so duplicated
#'   values are made unique via small decimal increments.
#' @export
check_for_slumps_add_decimals <- function(.ages) {
  # .ages <- df_chiro_codes[[8]]
  
  # check how many slumps
  .num_slumps <- length(unique(.ages$age_ce[duplicated(.ages$age_ce)]))
  
  if (.num_slumps != 0) {
    .id_slumps <- unique(.ages$age_ce[duplicated(.ages$age_ce)])
    if (unique(is.na(.id_slumps)) == FALSE) {
      for (.slump in 1:.num_slumps) {
        #     .slump <- 1
        .id_s <- .id_slumps[.slump]
        # replace
        .ages[.ages$age_ce == .id_s, "age_ce"] <-
          .ages[.ages$age_ce == .id_s, "age_ce"] +
          rev(seq(from = 0.01, by = 0.001, length.out = nrow(.ages[.ages$age_ce == .id_s, "age_ce"])))
      } # close loop
    }
  }
  
  res <- .ages
  
  return(res)
} # close function

#' Detect Number of Clusters from Broken-Stick Curve
#'
#' Uses the broken-stick criterion (dispersion vs bstick) to choose the
#' number of groups: returns the last \code{nGroups} before the first
#' crossing where dispersion drops below the broken-stick expectation.
#'
#' @param .bsticks A data frame with columns \code{dispersion}, \code{bstick}, and \code{nGroups}.
#' @return Integer: the suggested number of clusters (minimum 1).
#' @export
detect_num_clusts_bstick <- function(.bsticks) {
  # .bsticks <- all_bstick_chclusts[[5]]
  .get_sign <- .bsticks %>%
    mutate(sign = dispersion - bstick) %>%
    mutate(sign = sign(sign))
  .num_clusts <- .get_sign$nGroups[min(which(.get_sign$sign == -1)) - 1]
  
  if (identical(.num_clusts, integer(0)) == TRUE) {
    .num_clusts <- 1
  }
  
  return(.num_clusts)
}

#' Get Lake Names from Data Frames in a List
#'
#' Extracts the unique \code{lake} identifier from each data frame in a list.
#'
#' @param .lista A list of data frames, each with a \code{lake} column.
#' @return A list of vectors, one per element of \code{.lista}, giving unique lake names.
#' @export
get_names_list <- function(.lista) {
  .names <- lapply(.lista, function(x) unique(x[, "lake"]))
  return(.names)
}

#' Run Fuzzy C-Means (AMD) Clustering for Fixed Number of Clusters
#'
#' Runs fuzzy c-means clustering a set number of times and returns the
#' cluster assignment that maximizes mean adjusted membership (AMD).
#' Used when the number of clusters is already fixed (e.g. from \code{optAMDclusters}).
#'
#' @param .data Numeric matrix or data frame of features to cluster.
#' @param .iterations Number of clustering runs. Default 100.
#' @param .min_clusts Unused; kept for compatibility.
#' @param .opt_num_clusts Integer; number of clusters (required).
#' @return A one-column data frame of cluster labels (one per row of \code{.data}).
#' @importFrom e1071 cmeans
#' @export
getAMDclusters <- function(.data, .iterations = 100, .min_clusts = 2, .opt_num_clusts = NULL) {
  
  #.data <- select(list_master_table[[1]], -c(lake:core_depth_id))
  #.opt_num_clusts <- 6
  
  .maxpm <- c(0)
  .clst <- c()
  .probs <- c()
  .maxProbs <- c()
  
  for (k in 1:.iterations) { # loop over set number of iterations
    
    .clusts_k <- cmeans(.data, .opt_num_clusts, 20, verbose = F, method = "cmeans", m = 2) # fuzzy clustering
    
    .clusters <- as.numeric(.clusts_k$cluster) # assign each sample a cluster
    
    .prob <- .clusts_k$membership # get membership probabilities for each group
    
    .maxprob <- apply(.prob, 1, max) # max membership probability
    
    .mpm <- mean(.maxprob) - 1 / .opt_num_clusts # mean max membership probability
    
    # cat(i, "\n")
    
    if (.mpm > .maxpm) {
      .maxpm <- .mpm
      .clst <- .clusters
      .maxprobs <- .maxprob
      .probs <- .prob
      #print(.maxpm)
    }
  }
  
  .clst <- as.data.frame(.clst)
  
  # get n samples per cluster
  # table(.clst)
  
  return(.clst)
  
}

#' Select Optimal Number of Clusters Using AMD Heuristic
#'
#' This function estimates the optimal number of fuzzy clusters for a given dataset
#' using a heuristic based on Average Membership Degree (AMD). It iterates over a 
#' range of cluster sizes and repetitions, and selects the number of clusters that 
#' maximizes the mean adjusted maximum membership.
#'
#' @param .data A numeric matrix or data frame of features (e.g., species, traits) to be clustered.
#' @param .iterations Integer. Number of times to repeat the fuzzy clustering process. Default is 100.
#' @param .min_clusts Integer. Minimum number of clusters to test. Default is 2.
#' @param .num_groups Integer. Maximum number of clusters to test. Default is 12.
#'
#' @return An integer indicating the optimal number of clusters.
#'
#' @importFrom e1071 cmeans
#' @importFrom dplyr bind_cols
#' @export
optAMDclusters <- function(.data, .iterations = 100, .min_clusts = 2, .num_groups = 12) {
  
  .maxpms <- c()  # store permutation results across iterations
  
  for (k in 1:.iterations) {
    
    .maxpm <- c()  # store membership adjustment scores for each cluster size
    
    for (i in .min_clusts:.num_groups) {
      
      .clusts_k <- cmeans(.data, i, 20, verbose = FALSE, method = "cmeans", m = 2)
      .prob <- .clusts_k$membership  # membership probabilities
      
      # Extract the highest membership probability for each observation
      .maxprob <- apply(.prob, 1, max)
      
      # Compute adjusted mean membership (AMD heuristic)
      .mpm <- mean(.maxprob) - 1 / i
      .maxpm <- cbind(.maxpm, .mpm)
      
      n <- ncol(.maxpm)
      colnames(.maxpm)[n] <- paste("C", i, sep = "")  # label column with cluster count
    }
    
    .maxpms <- rbind(.maxpms, .maxpm)  # accumulate results for each iteration
  }
  
  # Helper function: maximum value by column
  colMax <- function(.col_data) {
    apply(.col_data, MARGIN = 2, max)
  }
  
  # Summary of results
  .maxpmean <- colMax(.maxpms)  # maximum AMD per cluster count
  .results <- bind_cols(.clusters = c(.min_clusts:.num_groups), .maxpmean = .maxpmean)
  
  .maxAMD <- round(max(.results$.maxpmean), 2)  # best AMD value
  .opt_clust_numb <- .results$.clusters[which(.results$.maxpmean == max(.results$.maxpmean))]  # optimal cluster count
  
  .output <- list(optimal = .opt_clust_numb)
  
  return(.output$optimal)
}

#' Plot AMD curve and extract optimal number of fuzzy clusters
#'
#' This function applies fuzzy clustering over a range of cluster numbers and multiple iterations,
#' computing Adjusted Membership Degree (AMD) scores to evaluate clustering quality.
#' It returns both the optimal number of clusters and a data frame summarizing AMD scores across cluster sizes.
#'
#' @param .data A numeric matrix or data frame containing the data to be clustered.
#' @param .iterations Integer. Number of iterations to perform the clustering process. Default is 100.
#' @param .min_clusts Integer. Minimum number of clusters to test. Default is 2.
#' @param .num_groups Integer. Maximum number of clusters to test. Default is 12.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{optimal}{The optimal number of clusters based on the maximum AMD score.}
#'   \item{results}{A data frame with the average AMD score for each number of clusters.}
#' }
#'
#' @seealso [e1071::cmeans()]
#' @export
PlotoptAMDclusters <- function(.data, .iterations = 100, .min_clusts = 2, .num_groups = 12) {
  
  .maxpms <- c()  # vector to store AMD scores across iterations
  
  for (k in 1:.iterations) {
    
    .maxpm <- c()  # vector to store AMD scores for each cluster count in this iteration
    
    for (i in .min_clusts:.num_groups) {
      
      .clusts_k <- cmeans(.data, i, 20, verbose = FALSE, method = "cmeans", m = 2)
      .prob <- .clusts_k$membership  # membership probabilities per cluster
      
      .maxprob <- apply(.prob, 1, max)  # highest membership per sample
      .mpm <- mean(.maxprob) - 1 / i  # Adjusted Membership Degree (AMD) heuristic
      
      .maxpm <- cbind(.maxpm, .mpm)  # store AMD for this cluster count
      n <- ncol(.maxpm)
      colnames(.maxpm)[n] <- paste("C", i, sep = "")
    }
    
    .maxpms <- rbind(.maxpms, .maxpm)  # accumulate results across iterations
  }
  
  # Helper function to get column-wise maximum
  colMax <- function(.col_data) {
    apply(.col_data, MARGIN = 2, max)
  }
  
  .maxpmean <- colMax(.maxpms)  # best AMD scores across iterations
  .results <- bind_cols(.clusters = c(.min_clusts:.num_groups), .maxpmean = .maxpmean)
  
  .maxAMD <- round(max(.results$.maxpmean), 2)
  .opt_clust_numb <- .results$.clusters[which(.results$.maxpmean == max(.results$.maxpmean))]
  
  .output <- list(
    optimal = .opt_clust_numb,
    results = .results
  )
  
  return(.output)
}

# ------------------------------------------------------------------------------
# Custom functions moved from main_script.Rmd
# ------------------------------------------------------------------------------

#' Fit GAM for Pollen Group with Site-Specific Smooths
#'
#' Fits a generalized additive model for percent abundance by age, with
#' optional site-specific smooths and random effect for site. Falls back
#' to a simple smooth if only one site is present.
#'
#' @param df Data frame with columns \code{percent}, \code{age_ce}, \code{sitename}, \code{group}.
#' @return A fitted \code{gam} object from \code{mgcv}.
#' @importFrom mgcv gam s
#' @export
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

#' Plot GAM Fits for a Pollen Group
#'
#' Builds a ggplot of observed percent vs age with fitted GAM curve and
#' annotation for adjusted R-squared and deviance explained.
#'
#' @param group_name Character; name of the group (e.g. "Tree", "Herb").
#' @param mod A fitted \code{gam} object from \code{fit_gam_group}.
#' @return A \code{ggplot} object.
#' @importFrom broom augment
#' @export
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

#' Map Band y-Coordinates to Original GAM (Vegetation) Scale
#'
#' Inverse transformation from discrete band y-positions to the vegetation
#' GAM scale, for use as a secondary axis in collage plots. Uses \code{ymin_band},
#' \code{ymax_band}, \code{veg_min}, \code{veg_max} from the calling environment.
#'
#' @param y Numeric; values in band space (e.g. 0.5, 1.5, 2.5, 3.5).
#' @return Numeric vector on vegetation scale.
#' @export
inv_gam_from_band <- function(y) {
  scales::rescale(y, from = c(ymin_band, ymax_band), to = c(veg_min, veg_max))
}

#' Bin Time Series into Fixed-Width Bins and Aggregate
#'
#' Bins a data frame by a time column (e.g. \code{age_ce}) using \code{bin_size}
#' from the calling environment, then summarizes specified columns by mean per bin.
#'
#' @param df Data frame with a time column.
#' @param time_col Name of the time column. Default \code{"age_ce"}.
#' @param cols Character vector of column names to average within each bin.
#' @return A tibble with \code{age_ce} (bin midpoints) and summarized \code{cols}.
#' @export
bin30 <- function(df, time_col = "age_ce", cols) {
  df %>%
    filter(!is.na(.data[[time_col]])) %>%
    mutate(age_bin = floor(.data[[time_col]] / bin_size) * bin_size) %>%
    group_by(age_bin) %>%
    summarise(across(all_of(cols), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    rename(age_ce = age_bin) %>%
    arrange(age_ce)
}

#' Compute Pure and Shared Variance (R²) for a Single Phase
#'
#' Runs redundancy analysis to partition variance in community data explained
#' by \code{selected_variables} into pure effects (each variable conditional on
#' the others) and shared variance. Uses \code{selected_variables} from the
#' calling environment.
#'
#' @param data Data frame with guild columns (\code{high_profile}:\code{predator})
#'   and \code{selected_variables} columns.
#' @param phase_name Character; label for the phase (e.g. "Phase 1").
#' @return A tibble with columns \code{phase}, \code{component}, \code{value}.
#' @importFrom vegan rda RsquareAdj
#' @export
calculate_pure_shared_phase <- function(data, phase_name) {
  comm <- data %>% dplyr::select(high_profile:predator) %>% dplyr::filter(complete.cases(.))
  env  <- data %>% dplyr::select(all_of(selected_variables)) %>% dplyr::filter(complete.cases(.))
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
  
  dplyr::bind_rows(pure_effects, shared)
}

#' Variance Partitioning by Phase (Pure A, Pure B, Shared)
#'
#' For each phase defined in \code{phase_breaks}, subsets the joined data and
#' runs RDA to compute pure effects of two predictors (e.g. NAO and vegetation)
#' and their shared variance. Uses \code{guild_cols}, \code{env_predictors},
#' \code{phase_breaks} from the calling environment by default.
#'
#' @param joined_df Data frame with \code{age_ce}, guild columns, and env predictor columns.
#' @param guild_cols Character vector of community column names.
#' @param env_predictors Length-two character vector (e.g. NAO and vegetation).
#' @param phase_breaks Named list of phases, each \code{c(lo, hi)} bounds for \code{age_ce}.
#' @param min_n Minimum rows per phase; otherwise returns NA for that phase.
#' @return A tibble with \code{phase}, \code{component}, \code{value}, \code{n_bins}.
#' @importFrom vegan rda RsquareAdj
#' @export
calc_phase_varpart <- function(joined_df,
                               guild_cols = guild_cols,
                               env_predictors = env_predictors,
                               phase_breaks = phase_breaks,
                               min_n = 5) {

  A <- env_predictors[1]
  B <- env_predictors[2]

  one_phase <- function(bounds, phase_name) {
    lo <- bounds[1]; hi <- bounds[2]

    sub <- joined_df %>%
      dplyr::filter(age_ce >= lo, age_ce < hi) %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(c(guild_cols, A, B)), ~ !is.na(.x)))

    if (nrow(sub) < min_n) {
      comps <- c(paste0("Pure ", A), paste0("Pure ", B), "Shared")
      return(tibble(
        phase     = phase_name,
        component = comps,
        value     = rep(NA_real_, length(comps)),
        n_bins    = rep(nrow(sub), length(comps))
      ))
    }

    comm <- sub %>% dplyr::select(dplyr::all_of(guild_cols)) %>% as.matrix()
    envA  <- sub %>% dplyr::select(dplyr::all_of(A)) %>% as.data.frame()
    envB  <- sub %>% dplyr::select(dplyr::all_of(B)) %>% as.data.frame()
    envAB <- sub %>% dplyr::select(dplyr::all_of(c(A, B))) %>% as.data.frame()

    mod_A <- vegan::rda(comm, envA, envB)
    pure_A <- max(0, vegan::RsquareAdj(mod_A)$adj.r.squared)

    mod_B <- vegan::rda(comm, envB, envA)
    pure_B <- max(0, vegan::RsquareAdj(mod_B)$adj.r.squared)

    mod_full <- vegan::rda(comm, envAB)
    total_r2 <- max(0, vegan::RsquareAdj(mod_full)$adj.r.squared)

    shared <- max(0, total_r2 - pure_A - pure_B)

    tibble(
      phase     = phase_name,
      component = c(paste0("Pure ", A), paste0("Pure ", B), "Shared"),
      value     = c(pure_A, pure_B, shared),
      n_bins    = nrow(sub)
    )
  }

  purrr::imap_dfr(phase_breaks, one_phase)
}

#' Build Phase Breaks List from Four Boundary Years
#'
#' Creates a named list of five phases with \code{age_ce} intervals:
#' Phase 1 = (-Inf, b1), Phase 2 = (b1, b2), ..., Phase 5 = (b4, Inf).
#'
#' @param b1,b2,b3,b4 Numeric; boundary years (e.g. 750, 1050, 1450, 1750).
#' @return Named list of five elements, each a length-two vector \code{c(lo, hi)}.
#' @export
make_phase_breaks <- function(b1, b2, b3, b4) {
  list(
    "Phase 1" = c(-Inf, b1),
    "Phase 2" = c(b1,  b2),
    "Phase 3" = c(b2,  b3),
    "Phase 4" = c(b3,  b4),
    "Phase 5" = c(b4,  Inf)
  )
}

#' Run Phase Variance Partition for One Boundary Set (Sensitivity)
#'
#' For a single quadruple of phase boundaries, runs \code{calc_phase_varpart}
#' and returns results in long form bound to grid metadata. Uses
#' \code{joined_df_30yr}, \code{guild_cols}, \code{env_predictors}, \code{min_n}
#' from the calling environment.
#'
#' @param b1,b2,b3,b4 Numeric; phase boundary years.
#' @param grid_id Identifier for this boundary set (e.g. row index).
#' @return A tibble with \code{grid_id}, boundaries, \code{any_phase_failed}, and long \code{phase}/\code{component}/\code{value}.
#' @export
run_one_boundary_set <- function(b1, b2, b3, b4, grid_id) {
  pb <- make_phase_breaks(b1, b2, b3, b4)

  out <- calc_phase_varpart(
    joined_df_30yr,
    guild_cols = guild_cols,
    env_predictors = env_predictors,
    phase_breaks = pb,
    min_n = min_n
  )

  # wide per phase
  out_wide <- out %>%
    dplyr::select(phase, component, value, n_bins) %>%
    tidyr::pivot_wider(names_from = component, values_from = value)

  # flag runs where any phase fell below min_n (will have NA values)
  any_phase_failed <- any(is.na(out_wide[[paste0("Pure ", env_predictors[1])]]) |
                            is.na(out_wide[[paste0("Pure ", env_predictors[2])]]))

  tibble(
    grid_id = grid_id,
    b1 = b1, b2 = b2, b3 = b3, b4 = b4,
    any_phase_failed = any_phase_failed
  ) %>%
    dplyr::bind_cols(
      out_wide %>%
        dplyr::arrange(phase) %>%
        # keep phase-level values in long form (more flexible)
        tidyr::pivot_longer(
          cols = -c(phase, n_bins),
          names_to = "component",
          values_to = "value"
        )
    )
}

#' Build Joined Regional Table from Biology and Environment Binned Data
#'
#' Left-joins regional biology (\code{bio_reg_30}) with binned NAO, NHST,
#' vegetation, and community proportion tables on \code{age_ce}, then
#' z-scores the selected environment columns.
#'
#' @param bio_reg_30 Data frame with \code{age_ce} and guild columns.
#' @param nao_30,nhst_30,veg_30,prop_comm_30 Data frames with \code{age_ce} and env/community columns.
#' @param env_cols Character vector of columns to scale. Default \code{c("NAO_Median_Value","Rmean","estimate")}.
#' @return A tibble with all joined columns and \code{env_cols} scaled.
#' @export
build_joined_from_bio_reg <- function(bio_reg_30,
                                      nao_30, nhst_30, veg_30, prop_comm_30,
                                      env_cols = c("NAO_Median_Value","Rmean","estimate")) {

  out <- bio_reg_30 %>%
    dplyr::left_join(nao_30, by = "age_ce") %>%
    dplyr::left_join(nhst_30, by = "age_ce") %>%
    dplyr::left_join(veg_30,  by = "age_ce") %>%
    dplyr::left_join(prop_comm_30, by = "age_ce")

  out %>%
    dplyr::mutate(dplyr::across(all_of(env_cols), ~ as.numeric(scale(.x))))
}

#' Moving-Window Variance Partition on Joined Regional Data
#'
#' For each sliding window of \code{window_size} years, computes RDA-based
#' variance partition (pure A, pure B, shared) for community data vs two
#' env predictors. Returns one row per window per component.
#'
#' @param joined_df_30yr Data frame with \code{age_ce}, \code{n_lakes}, guild columns, and env predictors.
#' @param window_size Window width in years. Default 300.
#' @param step_size Step between window starts. Default 30.
#' @param min_bins Minimum rows per window; otherwise that window is skipped.
#' @param guild_cols Character vector of community column names.
#' @param env_predictors Length-two character vector (e.g. \code{c("NAO_Median_Value", "estimate")}).
#' @return A tibble with \code{window_start}, \code{window_end}, \code{midpoint}, \code{component}, \code{value}, \code{n_lakes_mean}, \code{n_bins}.
#' @importFrom vegan rda RsquareAdj
#' @export
calc_window_varpart_joined <- function(joined_df_30yr,
                                       window_size = 300, step_size = 30,
                                       min_bins = 5,
                                       guild_cols = c("high_profile","low_profile","motile","euplanktonic",
                                                     "algivore","detritivore","plantivore","predator"),
                                       env_predictors = c("NAO_Median_Value", "estimate")) {

  stopifnot(length(env_predictors) == 2)
  A <- env_predictors[1]
  B <- env_predictors[2]

  dat <- joined_df_30yr %>%
    dplyr::select(age_ce, n_lakes, dplyr::all_of(guild_cols), dplyr::all_of(env_predictors)) %>%
    dplyr::arrange(age_ce)

  min_year <- floor(min(dat$age_ce, na.rm = TRUE))
  max_year <- ceiling(max(dat$age_ce, na.rm = TRUE))
  window_starts <- seq(min_year, max_year - window_size, by = step_size)

  one_window <- function(ws) {
    we <- ws + window_size
    sub <- dat %>% dplyr::filter(age_ce >= ws, age_ce < we)

    # require complete predictors within window
    sub <- sub %>% dplyr::filter(!is.na(.data[[A]]), !is.na(.data[[B]]))
    if (nrow(sub) < min_bins) return(NULL)

    # Response (community)
    comm <- sub %>%
      dplyr::select(dplyr::all_of(guild_cols)) %>%
      as.matrix()

    # Predictors
    envA <- sub %>% dplyr::select(dplyr::all_of(A)) %>% as.data.frame()
    envB <- sub %>% dplyr::select(dplyr::all_of(B)) %>% as.data.frame()
    envAB <- sub %>% dplyr::select(dplyr::all_of(c(A, B))) %>% as.data.frame()

    # Pure A (A | B): rda(X, Y, Z)
    mod_A <- vegan::rda(comm, envA, envB)
    pure_A <- max(0, vegan::RsquareAdj(mod_A)$adj.r.squared)

    # Pure B (B | A)
    mod_B <- vegan::rda(comm, envB, envA)
    pure_B <- max(0, vegan::RsquareAdj(mod_B)$adj.r.squared)

    # Total (A + B)
    mod_full <- vegan::rda(comm, envAB)
    total_r2 <- max(0, vegan::RsquareAdj(mod_full)$adj.r.squared)

    shared <- max(0, total_r2 - pure_A - pure_B)

    tibble(
      window_start = ws,
      window_end   = we,
      midpoint     = (ws + we) / 2,
      component    = c(paste0("Pure ", A), paste0("Pure ", B), "Shared"),
      value        = c(pure_A, pure_B, shared),
      n_lakes_mean = mean(sub$n_lakes, na.rm = TRUE),
      n_bins       = nrow(sub)
    )
  }

  purrr::map_dfr(window_starts, one_window)
}

#' Permutation p-Values for Moving-Window Variance Partition
#'
#' For one time window, runs full and partial RDA models and returns
#' permutation p-values for the full model, pure A, pure B, and sequential
#' terms. Uses \code{guild_cols} from the calling environment by default.
#'
#' @param data Data frame with \code{age_ce}, \code{n_lakes}, guild and env columns.
#' @param window_start Start year of the window.
#' @param window_size Window width. Default 300.
#' @param guild_cols Character vector of community column names.
#' @param A,B Character; names of the two env predictors. Defaults NAO and \code{estimate}.
#' @param min_n Minimum rows in window; otherwise returns NA p-values.
#' @param n_perm Number of permutations. Default 9999.
#' @return A one-row tibble with \code{window_start}, \code{window_end}, \code{midpoint}, \code{n_bins}, \code{mean_n_lakes}, and p-value columns.
#' @importFrom vegan rda
#' @export
calc_window_significance <- function(data, window_start,
                                     window_size = 300,
                                     guild_cols = guild_cols,
                                     A = "NAO_Median_Value",
                                     B = "estimate",
                                     min_n = 5,
                                     n_perm = 9999) {

  window_end <- window_start + window_size

  sub <- data %>%
    dplyr::select(age_ce, n_lakes, dplyr::all_of(guild_cols), dplyr::all_of(c(A, B))) %>%
    dplyr::filter(age_ce >= window_start, age_ce < window_end) %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(c(guild_cols, A, B)), ~ !is.na(.x)))

  n_bins <- nrow(sub)
  mean_n_lakes <- mean(sub$n_lakes, na.rm = TRUE)

  if (n_bins < min_n) {
    return(tibble(
      window_start = window_start,
      window_end   = window_end,
      midpoint     = (window_start + window_end) / 2,
      n_bins       = n_bins,
      mean_n_lakes = mean_n_lakes,
      p_full         = NA_real_,
      p_pure_A       = NA_real_,
      p_pure_B       = NA_real_,
      p_seq_A_then_B = NA_real_,
      p_seq_B_then_A = NA_real_
    ))
  }

  comm <- sub %>% dplyr::select(dplyr::all_of(guild_cols)) %>% as.matrix()
  env  <- sub %>% dplyr::select(dplyr::all_of(c(A, B))) %>% as.data.frame()

  # Full model: comm ~ A + B
  mod_full <- vegan::rda(as.formula(paste("comm ~", paste(c(A, B), collapse = " + "))), data = env)
  a_full <- safe_anova_p(mod_full, permutations = n_perm, return_table = TRUE)
  p_full <- if (is.null(a_full)) NA_real_ else as.numeric(a_full$`Pr(>F)`[1])

  # Pure A (A | B)
  mod_pure_A <- vegan::rda(as.formula(paste("comm ~", A, "+ Condition(", B, ")")), data = env)
  a_pA <- safe_anova_p(mod_pure_A, permutations = n_perm, return_table = TRUE)
  p_pure_A <- if (is.null(a_pA)) NA_real_ else as.numeric(a_pA$`Pr(>F)`[1])

  # Pure B (B | A)
  mod_pure_B <- vegan::rda(as.formula(paste("comm ~", B, "+ Condition(", A, ")")), data = env)
  a_pB <- safe_anova_p(mod_pure_B, permutations = n_perm, return_table = TRUE)
  p_pure_B <- if (is.null(a_pB)) NA_real_ else as.numeric(a_pB$`Pr(>F)`[1])

  # Sequential diagnostics: test SECOND term given FIRST (order-dependent)
  mod_A_then_B <- vegan::rda(as.formula(paste("comm ~", A, "+", B)), data = env)
  a_A_then_B <- safe_anova_p(mod_A_then_B, permutations = n_perm, by = "term", return_table = TRUE)
  p_seq_A_then_B <- if (is.null(a_A_then_B) || nrow(a_A_then_B) < 2) NA_real_ else as.numeric(a_A_then_B$`Pr(>F)`[2])

  mod_B_then_A <- vegan::rda(as.formula(paste("comm ~", B, "+", A)), data = env)
  a_B_then_A <- safe_anova_p(mod_B_then_A, permutations = n_perm, by = "term", return_table = TRUE)
  p_seq_B_then_A <- if (is.null(a_B_then_A) || nrow(a_B_then_A) < 2) NA_real_ else as.numeric(a_B_then_A$`Pr(>F)`[2])

  tibble(
    window_start = window_start,
    window_end   = window_end,
    midpoint     = (window_start + window_end) / 2,
    n_bins       = n_bins,
    mean_n_lakes = mean_n_lakes,
    p_full         = p_full,
    p_pure_A       = p_pure_A,
    p_pure_B       = p_pure_B,
    p_seq_A_then_B = p_seq_A_then_B,
    p_seq_B_then_A = p_seq_B_then_A
  )
}

#' Circularly Shift Age Bins Within Each Lake (Null Model)
#'
#' For each lake, randomly shifts the \code{age_ce} values so that temporal
#' order is preserved but alignment across lakes is broken. Used to generate
#' null distributions for variance partition statistics.
#'
#' @param bio_lake_30 Data frame with \code{lake} and \code{age_ce} (one row per bin per lake).
#' @return Same structure with \code{age_ce} circularly permuted within each \code{lake}.
#' @export
circular_shift_by_lake <- function(bio_lake_30) {
  bio_lake_30 %>%
    dplyr::group_by(lake) %>%
    dplyr::arrange(age_ce, .by_group = TRUE) %>%
    dplyr::group_modify(~{
      n <- nrow(.x)
      if (n <= 1) return(.x)
      shift <- sample.int(n, 1) - 1
      .x %>% dplyr::mutate(age_ce = age_ce[((dplyr::row_number() + shift - 1) %% n) + 1])
    }) %>%
    dplyr::ungroup()
}

#' Plot Observed vs Null Envelope for One Variance Component
#'
#' Plots observed adjusted R² over time (midpoint) for a given component
#' together with the null median (dashed). Uses \code{r2_base} and
#' \code{null_env} from the calling environment.
#'
#' @param component_name Character; component to plot (e.g. \code{"Pure estimate"}, \code{"Pure NAO_Median_Value"}).
#' @return A \code{ggplot} object.
#' @export
plot_obs_vs_null <- function(component_name) {

  obs <- r2_base %>%
    dplyr::filter(component == component_name) %>%
    dplyr::select(midpoint, obs = value)

  nul <- null_env %>%
    dplyr::filter(component == component_name) %>%
    dplyr::select(midpoint, q05, q50, q95)

  df <- obs %>% dplyr::inner_join(nul, by = "midpoint")

  ggplot(df, aes(x = midpoint)) +
    geom_line(aes(y = q50),
              linewidth = 0.8, linetype = "dashed", color = "grey30") +
    geom_line(aes(y = obs),
              linewidth = 1.2, color = "black") +
    labs(
      x = "Year (CE)",
      y = "Adjusted R²",
      title = paste0(component_name, ": observed vs null")
    ) +
    theme_minimal(base_size = 13)
}

#' Plot Observed vs Null with Points Where Observed Exceeds Null
#'
#' Same as \code{plot_obs_vs_null} but with optional custom \code{ylab} and
#' \code{title_text}, and points highlighted in red where observed > null median.
#' Uses \code{r2_base} and \code{null_env} from the calling environment.
#'
#' @param component_name Character; component to plot.
#' @param ylab Character; y-axis label.
#' @param title_text Character; plot title.
#' @return A \code{ggplot} object.
#' @export
plot_obs_vs_null_points <- function(component_name, ylab, title_text) {

  obs <- r2_base %>%
    dplyr::filter(component == component_name) %>%
    dplyr::select(midpoint, obs = value)

  nul <- null_env %>%
    dplyr::filter(component == component_name) %>%
    dplyr::select(midpoint, q50)

  df <- obs %>%
    dplyr::inner_join(nul, by = "midpoint") %>%
    dplyr::mutate(above_null = obs > q50)

  ggplot(df, aes(x = midpoint)) +
    # Null median (dashed)
    geom_line(
      aes(y = q50),
      linewidth = 0.8,
      linetype = "dashed",
      color = "grey40"
    ) +
    # Observed (black line)
    geom_line(
      aes(y = obs),
      linewidth = 1.1,
      color = "black"
    ) +
    # Points where observed > null median
    geom_point(
      data = df %>% dplyr::filter(above_null),
      aes(y = obs),
      color = "red",
      size = 2
    ) +
    labs(
      x = "Year (CE)",
      y = ylab,
      title = title_text
    ) +
    theme_minimal(base_size = 13)
}

#' Map Vegetation Scale to DCA Axis Scale
#'
#' Linear rescaling from a fixed vegetation range to the DCA axis range.
#' Uses \code{veg_limits} and \code{dca_rng} from the calling environment.
#'
#' @param v Numeric; value(s) on vegetation scale.
#' @return Numeric; value(s) on DCA scale.
#' @export
veg_to_dca <- function(v) {
  (v - veg_limits[1]) / diff(veg_limits) * diff(dca_rng) + dca_rng[1]
}

#' Map DCA Axis Scale to Vegetation Scale
#'
#' Inverse of \code{veg_to_dca}: rescales from DCA range to vegetation range.
#' Uses \code{dca_rng} and \code{veg_limits} from the calling environment.
#'
#' @param y Numeric; value(s) on DCA scale.
#' @return Numeric; value(s) on vegetation scale.
#' @export
dca_to_veg <- function(y) {
  (y - dca_rng[1]) / diff(dca_rng) * diff(veg_limits) + veg_limits[1]
}

#' Compute Pure and Shared Variance (R²) for One Moving Window
#'
#' For a single time window, partitions variance in community data explained
#' by \code{selected_variables} into pure effects and shared, via RDA.
#' Returns \code{NULL} if fewer than 5 complete rows in the window.
#'
#' @param data Data frame with guild columns and \code{selected_variables} columns.
#' @param window_start Start year of the window.
#' @param window_size Window width in years.
#' @param selected_variables Character vector of env variable names (e.g. \code{c("NAO_Median_Value", "estimate")}).
#' @return A tibble with \code{window_start}, \code{window_end}, \code{component}, \code{value}, or \code{NULL}.
#' @importFrom vegan rda RsquareAdj
#' @export
calculate_pure_shared_window <- function(data, window_start, window_size, selected_variables) {
  window_end <- window_start + window_size
  subset <- data %>% dplyr::filter(age_ce >= window_start, age_ce < window_end)

  comm <- subset %>% dplyr::select(high_profile:predator) %>% dplyr::filter(complete.cases(.))
  env  <- subset %>% dplyr::select(all_of(selected_variables)) %>% dplyr::filter(complete.cases(.))
  n <- min(nrow(comm), nrow(env))

  if (n < 5) return(NULL)

  comm <- comm[1:n, ]
  env  <- env[1:n, ]

  pure_effects <- purrr::map_dfr(selected_variables, function(var) {
    others <- setdiff(selected_variables, var)
    formula <- as.formula(
      paste0("comm ~ ", var, " + Condition(", paste(others, collapse = " + "), ")")
    )
    mod <- vegan::rda(formula = formula, data = env)
    tibble(
      window_start = window_start,
      window_end   = window_end,
      component    = ifelse(var == "estimate", "Pure Vegetation", "Pure Climate"),
      value        = max(0, vegan::RsquareAdj(mod)$adj.r.squared)
    )
  })

  mod_full <- vegan::rda(comm ~ ., data = env)
  total_r2 <- max(0, vegan::RsquareAdj(mod_full)$adj.r.squared)
  shared_val <- total_r2 - sum(pure_effects$value)

  shared <- tibble(
    window_start = window_start,
    window_end   = window_end,
    component    = "Shared",
    value        = max(0, shared_val)
  )

  dplyr::bind_rows(pure_effects, shared)
}

#' Safe Permutation ANOVA p-Value or Table
#'
#' Wraps \code{anova(..., permutations)} in \code{tryCatch}. Returns a single
#' p-value (for \code{which_row}), the full ANOVA table, or \code{NA}/\code{NULL}
#' on failure. Used with \code{vegan} ordination models.
#'
#' @param mod A model object with an \code{anova} method (e.g. \code{rda}).
#' @param permutations Number of permutations. Default 9999.
#' @param by Passed to \code{anova} (e.g. \code{"term"} for sequential tests).
#' @param which_row Row index of the ANOVA table for the p-value when \code{return_table} is FALSE. Default 1.
#' @param return_table If \code{TRUE}, return the full ANOVA table; otherwise return a single numeric p-value.
#' @return Numeric p-value, data frame (ANOVA table), or \code{NA_real_}/\code{NULL} on error.
#' @export
safe_anova_p <- function(mod, permutations = 9999, by = NULL, which_row = 1, return_table = FALSE) {
  out <- tryCatch(stats::anova(mod, permutations = permutations, by = by), error = function(e) NULL)
  if (is.null(out)) return(if (return_table) NULL else NA_real_)
  if (return_table) return(out)
  if (!"Pr(>F)" %in% colnames(out)) return(NA_real_)
  if (nrow(out) < which_row) return(NA_real_)
  as.numeric(out$`Pr(>F)`[which_row])
}

#' Prepare long-format data for a lake stratiplot (top-N taxa + Other)
#'
#' Aggregates duplicate \code{age_ce} levels by summing counts, drops all-zero
#' taxa, keeps the \code{n_top} most abundant species, and returns relative
#' abundances interpolated onto a fine age grid for continuous stacked bands.
#'
#' @param abund_wide Data frame of species counts with optional metadata columns.
#' @param codes_df Sample metadata including \code{age_ce}, aligned row-wise with
#'   \code{abund_wide}.
#' @param meta_cols Optional extra column names to exclude from the species matrix.
#'   Common metadata columns (\code{lake}, \code{core_depth_id}, \code{depth_top},
#'   \code{dec_depth}, \code{age_ce}, etc.) are always dropped automatically.
#' @param n_top Number of taxa to retain before pooling the remainder as \code{Other}.
#' @param age_step Age increment (years CE) for the interpolated grid. When
#'   \code{NULL}, uses one quarter of the median sample spacing (minimum 1 year).
#' @return A list with \code{data} (long tibble), \code{top_taxa}, and
#'   \code{age_range}.
#' @export
prepare_lake_stratiplot_data <- function(
    abund_wide,
    codes_df,
    meta_cols = NULL,
    n_top = 10L,
    age_step = NULL) {
  default_meta_cols <- c(
    "lake", "core_depth_id", "core_id", "sec_cor", "samp_id", "samp_dep",
    "dec_depth", "depth_top", "depth_bot", "age_bp", "age_ce"
  )
  drop_cols <- unique(c(default_meta_cols, meta_cols))

  abund <- abund_wide %>%
    dplyr::select(-dplyr::any_of(drop_cols)) %>%
    dplyr::mutate(.row_id = dplyr::row_number())

  meta <- codes_df %>%
    dplyr::mutate(.row_id = dplyr::row_number()) %>%
    dplyr::select(.row_id, age_ce) %>%
    dplyr::filter(!is.na(age_ce))

  combined <- abund %>%
    dplyr::inner_join(meta, by = ".row_id") %>%
    dplyr::select(-.row_id)

  combined <- combined %>%
    dplyr::group_by(age_ce) %>%
    dplyr::summarise(
      dplyr::across(dplyr::where(is.numeric), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::arrange(age_ce)

  species_cols <- setdiff(names(combined), "age_ce")
  species_cols <- species_cols[colSums(combined[species_cols], na.rm = TRUE) > 0]

  abund_mat <- as.matrix(combined[, species_cols, drop = FALSE])
  totals <- colSums(abund_mat, na.rm = TRUE)
  top_taxa <- names(sort(totals, decreasing = TRUE))[seq_len(min(n_top, length(totals)))]
  other_cols <- setdiff(species_cols, top_taxa)

  row_tot <- rowSums(abund_mat, na.rm = TRUE)
  row_tot[row_tot == 0] <- NA_real_

  pct_top <- sweep(abund_mat[, top_taxa, drop = FALSE], 1, row_tot, "/") * 100
  pct_other <- if (length(other_cols) > 0) {
    rowSums(abund_mat[, other_cols, drop = FALSE], na.rm = TRUE) / row_tot * 100
  } else {
    rep(0, length(row_tot))
  }

  plot_wide <- dplyr::bind_cols(
    combined[, "age_ce", drop = FALSE],
    as.data.frame(pct_top),
    Other = pct_other
  )

  taxon_levels <- c(top_taxa, "Other")
  ages <- plot_wide$age_ce

  if (is.null(age_step)) {
    age_step <- if (length(ages) > 1L) {
      max(1, stats::median(diff(ages), na.rm = TRUE) / 4)
    } else {
      1
    }
  }

  fine_ages <- seq(min(ages), max(ages), by = age_step)
  interp_mat <- vapply(
    taxon_levels,
    function(tax) {
      stats::approx(
        x = ages,
        y = plot_wide[[tax]],
        xout = fine_ages,
        rule = 2
      )$y
    },
    FUN.VALUE = numeric(length(fine_ages))
  )
  interp_mat <- pmax(interp_mat, 0)
  row_sums <- rowSums(interp_mat)
  row_sums[row_sums == 0] <- 1
  interp_mat <- sweep(interp_mat, 1, row_sums, "/") * 100

  plot_long <- dplyr::as_tibble(interp_mat) %>%
    dplyr::mutate(age_ce = fine_ages) %>%
    tidyr::pivot_longer(-age_ce, names_to = "taxon", values_to = "percent") %>%
    dplyr::filter(!is.na(percent))

  plot_long$taxon <- factor(plot_long$taxon, levels = taxon_levels)

  list(
    data = plot_long,
    top_taxa = top_taxa,
    age_range = range(fine_ages, na.rm = TRUE)
  )
}

#' Continuous stacked stratiplot for one lake (age CE on y-axis; present at top)
#'
#' @param prep Output of \code{\link{prepare_lake_stratiplot_data}}.
#' @param lake_name Lake title.
#' @param taxon_group Subtitle (e.g. \code{"Diatoms"} or \code{"Chironomids"}).
#' @param base_size Base ggplot theme size.
#' @param show_y_axis If \code{FALSE}, hide the y-axis title/labels (for multi-panel layouts).
#' @param age_limits Optional numeric length-2 vector of shared y-axis limits.
#' @param panel_tag Optional panel tag (e.g. \code{"(a)"}).
#' @return A \code{ggplot} object.
#' @export
plot_lake_stratiplot <- function(
    prep,
    lake_name,
    taxon_group,
    base_size = 10,
    show_y_axis = TRUE,
    age_limits = NULL,
    panel_tag = NULL) {
  fill_vals <- c(
    stats::setNames(
      grDevices::hcl.colors(length(prep$top_taxa), palette = "Spectral"),
      prep$top_taxa
    ),
    Other = "grey85"
  )

  plot_data <- prep$data %>%
    dplyr::arrange(age_ce, taxon)

  y_lims <- if (is.null(age_limits)) prep$age_range else age_limits

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = percent, y = age_ce, fill = taxon)) +
    ggplot2::geom_area(
      position = "stack",
      orientation = "y",
      linewidth = 0
    ) +
    ggplot2::scale_fill_manual(values = fill_vals, name = NULL) +
    ggplot2::coord_cartesian(xlim = c(0, 100), ylim = y_lims, clip = "off") +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(
      tag = panel_tag,
      title = taxon_group,
      subtitle = "Top taxa (%)",
      x = "Relative abundance (%)",
      y = if (show_y_axis) "Age (CE)" else NULL
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.tag = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.subtitle = ggplot2::element_text(size = base_size - 1),
      panel.grid.minor = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 6),
      legend.key.size = grid::unit(0.28, "cm"),
      legend.position = "bottom"
    )

  if (!show_y_axis) {
    p <- p + ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )
  }

  p
}

#' Lake-level DCA1 GAM panel aligned with stratiplots (age on y; present at top)
#'
#' @param dca_df Data frame with \code{age_ce} and \code{DCA1}.
#' @param taxon_group Panel title (e.g. \code{"Diatoms"}).
#' @param line_colour Colour for points and smooth.
#' @param base_size Base ggplot theme size.
#' @param show_y_axis If \code{FALSE}, hide the y-axis title/labels.
#' @param age_limits Optional shared y-axis limits.
#' @param panel_tag Optional panel tag.
#' @return A \code{ggplot} object.
#' @export
plot_lake_dca_panel <- function(
    dca_df,
    taxon_group,
    line_colour = "#440154",
    base_size = 10,
    show_y_axis = TRUE,
    age_limits = NULL,
    panel_tag = NULL) {
  dca_df <- dca_df %>%
    dplyr::filter(!is.na(age_ce), !is.na(DCA1)) %>%
    dplyr::arrange(age_ce)

  y_lims <- if (is.null(age_limits)) {
    range(dca_df$age_ce, na.rm = TRUE)
  } else {
    age_limits
  }

  p <- ggplot2::ggplot(dca_df, ggplot2::aes(x = DCA1, y = age_ce)) +
    ggplot2::geom_point(alpha = 0.45, size = 1.4, colour = line_colour) +
    ggplot2::geom_smooth(
      method = "gam",
      formula = y ~ s(x, k = 8),
      se = TRUE,
      colour = line_colour,
      fill = line_colour,
      alpha = 0.2,
      linewidth = 1.1,
      orientation = "y"
    ) +
    ggplot2::coord_cartesian(ylim = y_lims) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.02))) +
    ggplot2::labs(
      tag = panel_tag,
      title = taxon_group,
      subtitle = "DCA1 (GAM)",
      x = "DCA1",
      y = if (show_y_axis) "Age (CE)" else NULL
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.tag = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.subtitle = ggplot2::element_text(size = base_size - 1),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (!show_y_axis) {
    p <- p + ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )
  }

  p
}

#' Four-panel lake figure: diatom/chironomid stratiplots + DCA GAMs
#'
#' Layout: diatom stratiplot | diatom DCA | chironomid stratiplot | chironomid DCA.
#'
#' @param lake Lake name (used internally).
#' @param diat_prep Prepared diatom stratiplot data.
#' @param chiro_prep Prepared chironomid stratiplot data.
#' @param diat_dca Diatom DCA data for the lake.
#' @param chiro_dca Chironomid DCA data for the lake.
#' @param lake_title Optional display title (e.g. \code{"Caldeirão [CR]"}).
#' @param base_size Base ggplot theme size.
#' @return A patchwork object.
#' @export
plot_lake_strati_dca_composite <- function(
    lake,
    diat_prep,
    chiro_prep,
    diat_dca,
    chiro_dca,
    lake_title = NULL,
    base_size = 9) {
  display_title <- if (is.null(lake_title) || !nzchar(lake_title)) lake else lake_title

  age_limits <- range(
    c(diat_prep$age_range, chiro_prep$age_range,
      diat_dca$age_ce, chiro_dca$age_ce),
    na.rm = TRUE
  )

  diat_col <- viridisLite::viridis(9)[9]
  chiro_col <- viridisLite::viridis(9)[1]

  p_diat_st <- plot_lake_stratiplot(
    diat_prep,
    lake_name = lake,
    taxon_group = "Diatoms",
    base_size = base_size,
    show_y_axis = TRUE,
    age_limits = age_limits,
    panel_tag = "(a)"
  )
  p_diat_dca <- plot_lake_dca_panel(
    diat_dca,
    taxon_group = "Diatoms",
    line_colour = diat_col,
    base_size = base_size,
    show_y_axis = FALSE,
    age_limits = age_limits,
    panel_tag = "(b)"
  )
  p_chiro_st <- plot_lake_stratiplot(
    chiro_prep,
    lake_name = lake,
    taxon_group = "Chironomids",
    base_size = base_size,
    show_y_axis = FALSE,
    age_limits = age_limits,
    panel_tag = "(c)"
  )
  p_chiro_dca <- plot_lake_dca_panel(
    chiro_dca,
    taxon_group = "Chironomids",
    line_colour = chiro_col,
    base_size = base_size,
    show_y_axis = FALSE,
    age_limits = age_limits,
    panel_tag = "(d)"
  )

  (p_diat_st | p_diat_dca | p_chiro_st | p_chiro_dca) +
    patchwork::plot_annotation(
      title = display_title,
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 14))
    ) +
    patchwork::plot_layout(widths = c(1.35, 0.85, 1.35, 0.85))
}

#' Export Appendix S2 Word document of lake stratiplot + DCA figures
#'
#' One landscape page per lake, preceded by a short explanatory cover page.
#'
#' @param plots Named list of ggplot/patchwork objects (names = lake keys).
#' @param lake_titles Named character vector mapping lake keys to display titles.
#' @param outfile Output \code{.docx} path.
#' @param n_top Number of top taxa shown in stratiplots (for the caption text).
#' @param fig_width Figure width in inches when rasterising.
#' @param fig_height Figure height in inches when rasterising.
#' @param dpi Raster resolution.
#' @return Invisibly, the path to \code{outfile}.
#' @export
export_appendix_s2_strati_dca <- function(
    plots,
    lake_titles = NULL,
    outfile = "supplementary/Appendix_S2.docx",
    n_top = 10L,
    fig_width = 14,
    fig_height = 7,
    dpi = 300) {
  if (!requireNamespace("officer", quietly = TRUE)) {
    stop("Package 'officer' is required to write Appendix S2.", call. = FALSE)
  }

  lake_keys <- names(plots)
  if (is.null(lake_keys) || any(!nzchar(lake_keys))) {
    stop("`plots` must be a named list of lake figures.", call. = FALSE)
  }

  title_of <- function(lake) {
    if (!is.null(lake_titles) && lake %in% names(lake_titles)) {
      unname(lake_titles[[lake]])
    } else {
      lake
    }
  }

  dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
  tmp_dir <- tempfile("appendix_s2_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  png_paths <- character(length(lake_keys))
  names(png_paths) <- lake_keys
  for (lake in lake_keys) {
    png_paths[[lake]] <- file.path(tmp_dir, paste0(gsub(" ", "_", lake), ".png"))
    ggplot2::ggsave(
      png_paths[[lake]],
      plot = plots[[lake]],
      width = fig_width,
      height = fig_height,
      dpi = dpi,
      bg = "white"
    )
  }

  intro <- paste0(
    "Appendix S2. Lake-level stratigraphic diagrams and community turnover ",
    "(diatoms and chironomids).\n\n",
    "Each subsequent page shows one study lake. Panels are: (a) diatom relative ",
    "abundances for the top ", n_top, " taxa (remaining taxa pooled as Other); ",
    "(b) diatom DCA1 generalised additive model (GAM) smooth with 95% confidence ",
    "band; (c) chironomid relative abundances for the top ", n_top, " taxa ",
    "(plus Other); (d) chironomid DCA1 GAM smooth with 95% confidence band. ",
    "Age (Common Era) is shared on the vertical axis, with the present at the top. ",
    "Lake titles use the same island abbreviations as Fig. 2b ",
    "(CR = Corvo, FL = Flores, PI = Pico, SM = São Miguel, TE = Terceira)."
  )

  doc <- officer::read_docx()
  doc <- officer::body_add_par(doc, "Appendix S2", style = "heading 1")
  for (para in strsplit(intro, "\n\n", fixed = TRUE)[[1]]) {
    doc <- officer::body_add_par(doc, para, style = "Normal")
  }

  # Landscape pages for the lake figures
  landscape_section <- officer::prop_section(
    page_size = officer::page_size(orient = "landscape", width = 11.69, height = 8.27),
    page_margins = officer::page_mar(
      bottom = 0.5, top = 0.5, right = 0.5, left = 0.5,
      header = 0.3, footer = 0.3, gutter = 0
    )
  )

  for (i in seq_along(lake_keys)) {
    lake <- lake_keys[[i]]
    doc <- officer::body_add_break(doc)
    doc <- officer::body_add_par(doc, title_of(lake), style = "heading 2")
    doc <- officer::body_add_img(
      doc,
      src = png_paths[[lake]],
      width = 10.6,
      height = 5.3
    )
  }

  # Apply landscape to the whole document after the first page break content
  doc <- officer::body_set_default_section(doc, landscape_section)

  print(doc, target = outfile)
  invisible(outfile)
}

#' Compute simultaneous CI derivatives for the regional age smooth
#'
#' @param mod Fitted HGAM containing \code{s(age_ce)}.
#' @param series_label Label for the series (e.g. \code{"Producers"}).
#' @param level Confidence level for simultaneous intervals. Default 0.9.
#' @param eps Finite-difference step for \code{gratia::derivatives}. Default 1 year.
#' @param seed RNG seed for simultaneous intervals.
#' @return Tibble with derivative estimates, CIs, and significance flags.
#' @export
compute_regional_gam_derivatives <- function(
    mod,
    series_label,
    level = 0.9,
    eps = 1,
    seed = 22) {
  dr <- gratia::derivatives(
    mod,
    select = "s(age_ce)",
    eps = eps,
    level = level,
    interval = "simultaneous",
    seed = seed
  )

  age_col <- if ("age_ce" %in% names(dr)) "age_ce" else {
    setdiff(names(dr), c(
      ".smooth", ".by", ".fs", ".derivative", ".se", ".crit",
      ".lower_ci", ".upper_ci"
    ))[1]
  }

  tibble::tibble(
    series = series_label,
    age_ce = as.numeric(dr[[age_col]]),
    derivative = as.numeric(dr$.derivative),
    lower_ci = as.numeric(dr$.lower_ci),
    upper_ci = as.numeric(dr$.upper_ci),
    significant = (dr$.lower_ci > 0) | (dr$.upper_ci < 0)
  ) %>%
    dplyr::arrange(age_ce)
}

#' Summarise significant-change periods from GAM derivatives
#'
#' Identifies contiguous intervals where the simultaneous CI for the derivative
#' excludes zero, then reports the first onset and the age of maximum absolute
#' derivative within significant periods.
#'
#' @param deriv_df Output of \code{\link{compute_regional_gam_derivatives}}.
#' @param min_run Minimum consecutive evaluation points to retain a period.
#' @return List with \code{periods}, \code{summary}, and \code{deriv}.
#' @export
summarise_gam_change_periods <- function(deriv_df, min_run = 3L) {
  d <- deriv_df %>%
    dplyr::arrange(age_ce) %>%
    dplyr::mutate(
      run_id = cumsum(significant != dplyr::lag(significant, default = FALSE))
    )

  periods <- d %>%
    dplyr::filter(significant, is.finite(age_ce)) %>%
    dplyr::group_by(series, run_id) %>%
    dplyr::reframe(
      start_ce = min(age_ce, na.rm = TRUE),
      end_ce = max(age_ce, na.rm = TRUE),
      n_points = dplyr::n(),
      mean_abs_deriv = mean(abs(derivative), na.rm = TRUE),
      peak_abs_deriv = max(abs(derivative), na.rm = TRUE),
      peak_ce = age_ce[which.max(abs(derivative))][1],
      direction = dplyr::if_else(
        mean(derivative, na.rm = TRUE) >= 0,
        "increase",
        "decrease"
      )
    ) %>%
    dplyr::filter(is.finite(start_ce), is.finite(end_ce), n_points >= min_run) %>%
    dplyr::arrange(start_ce) %>%
    dplyr::mutate(period_id = dplyr::row_number())

  summary_tbl <- if (nrow(periods) == 0) {
    tibble::tibble(
      series = unique(deriv_df$series),
      first_onset_ce = NA_real_,
      first_onset_after_800_ce = NA_real_,
      first_end_after_800_ce = NA_real_,
      first_direction = NA_character_,
      first_direction_after_800 = NA_character_,
      first_peak_ce = NA_real_,
      n_significant_periods = 0L,
      earliest_significant_ce = NA_real_,
      latest_significant_ce = NA_real_
    )
  } else {
    periods %>%
      dplyr::group_by(series) %>%
      dplyr::summarise(
        first_onset_ce = min(start_ce),
        first_onset_after_800_ce = {
          post_idx <- which(start_ce >= 800)
          if (length(post_idx) == 0) NA_real_ else min(start_ce[post_idx])
        },
        first_end_after_800_ce = {
          post_idx <- which(start_ce >= 800)
          if (length(post_idx) == 0) NA_real_ else end_ce[post_idx][which.min(start_ce[post_idx])]
        },
        first_direction = direction[which.min(start_ce)],
        first_direction_after_800 = {
          post_idx <- which(start_ce >= 800)
          if (length(post_idx) == 0) {
            NA_character_
          } else {
            direction[post_idx][which.min(start_ce[post_idx])]
          }
        },
        first_peak_ce = peak_ce[which.min(start_ce)],
        n_significant_periods = dplyr::n(),
        earliest_significant_ce = min(start_ce),
        latest_significant_ce = max(end_ce),
        .groups = "drop"
      )
  }

  list(periods = periods, summary = summary_tbl, deriv = d)
}

#' Fit a lake-level DCA1 ~ s(age_ce) GAM
#'
#' @param dca_df Data frame with \code{age_ce} and \code{DCA1}.
#' @param min_n Minimum number of finite observations required.
#' @param k Basis dimension for the age smooth (capped by sample size).
#' @return Fitted \code{gam} object, or \code{NULL} if fitting fails.
#' @export
fit_lake_dca_gam <- function(dca_df, min_n = 8L, k = 8L) {
  df <- dca_df %>%
    dplyr::transmute(
      age_ce = as.numeric(age_ce),
      DCA1 = as.numeric(DCA1)
    ) %>%
    dplyr::filter(is.finite(age_ce), is.finite(DCA1)) %>%
    dplyr::arrange(age_ce)

  if (nrow(df) < min_n) {
    return(NULL)
  }

  k_use <- min(as.integer(k), max(3L, floor(nrow(df) / 2L) - 1L))

  tryCatch(
    mgcv::gam(
      DCA1 ~ s(age_ce, bs = "tp", k = k_use),
      data = df,
      method = "REML"
    ),
    error = function(e) NULL
  )
}

#' Derivative timing for one lake x trophic group
#'
#' @param dca_df Lake-level DCA data.
#' @param lake_name Lake name.
#' @param group_label \code{"Producers"} or \code{"Consumers"}.
#' @param min_n Minimum observations to fit the GAM.
#' @param onset_min_ce Ignore significant onsets before this age.
#' @return List with \code{summary}, \code{periods}, \code{smooth}, and \code{deriv}.
#' @export
analyse_lake_gam_timing <- function(
    dca_df,
    lake_name,
    group_label,
    min_n = 8L,
    onset_min_ce = 800) {
  lake_df <- dca_df %>%
    dplyr::filter(.data$lake == .env$lake_name)

  n_obs <- nrow(lake_df %>% dplyr::filter(is.finite(age_ce), is.finite(DCA1)))
  mod <- fit_lake_dca_gam(lake_df, min_n = min_n)
  if (is.null(mod)) {
    empty <- tibble::tibble(
      lake = lake_name,
      group = group_label,
      n = n_obs,
      first_onset_ce = NA_real_,
      first_onset_after_min_ce = NA_real_,
      first_end_after_min_ce = NA_real_,
      first_direction = NA_character_,
      first_direction_after_min = NA_character_,
      first_peak_ce = NA_real_,
      total_dca_change = NA_real_,
      n_significant_periods = 0L
    )
    return(list(
      summary = empty,
      periods = tibble::tibble(),
      smooth = tibble::tibble(),
      deriv = tibble::tibble()
    ))
  }

  smooth <- gratia::smooth_estimates(mod) %>%
    dplyr::filter(.smooth == "s(age_ce)") %>%
    gratia::add_confint() %>%
    dplyr::transmute(
      lake = lake_name,
      group = group_label,
      age_ce = as.numeric(age_ce),
      estimate = as.numeric(.estimate),
      lower_ci = as.numeric(.lower_ci),
      upper_ci = as.numeric(.upper_ci)
    )

  deriv <- compute_regional_gam_derivatives(mod, group_label) %>%
    dplyr::mutate(lake = lake_name, group = group_label)

  change <- summarise_gam_change_periods(deriv)
  total_change <- diff(range(smooth$estimate, na.rm = TRUE))

  post <- change$periods %>%
    dplyr::filter(start_ce >= onset_min_ce)

  summary_tbl <- change$summary %>%
    dplyr::transmute(
      lake = lake_name,
      group = group_label,
      n = n_obs,
      first_onset_ce = first_onset_ce,
      first_onset_after_min_ce = if (nrow(post) == 0) NA_real_ else min(post$start_ce),
      first_end_after_min_ce = if (nrow(post) == 0) {
        NA_real_
      } else {
        post$end_ce[which.min(post$start_ce)]
      },
      first_direction = first_direction,
      first_direction_after_min = if (nrow(post) == 0) {
        NA_character_
      } else {
        post$direction[which.min(post$start_ce)]
      },
      first_peak_ce = first_peak_ce,
      total_dca_change = total_change,
      n_significant_periods = n_significant_periods
    )

  periods <- change$periods %>%
    dplyr::mutate(lake = lake_name, group = group_label)

  list(
    summary = summary_tbl,
    periods = periods,
    smooth = smooth,
    deriv = deriv
  )
}

#' Lake-level producer vs consumer GAM change timing
#'
#' @param diat_df Diatom DCA data with \code{lake}, \code{age_ce}, \code{DCA1}.
#' @param chiro_df Chironomid DCA data with \code{lake}, \code{age_ce}, \code{DCA1}.
#' @param lakes Optional vector of lakes to include.
#' @param exclude_lakes Lakes to drop.
#' @param min_n Minimum observations per lake/group.
#' @param onset_min_ce Minimum age for onset comparison.
#' @return List with \code{summary}, \code{comparison}, \code{periods}, \code{smooths}, \code{derivs}.
#' @export
compute_lake_level_gam_timing <- function(
    diat_df,
    chiro_df,
    lakes = NULL,
    exclude_lakes = c("Fogo", "Furnas"),
    min_n = 8L,
    onset_min_ce = 800) {
  lake_pool <- intersect(
    unique(diat_df$lake),
    unique(chiro_df$lake)
  )
  lake_pool <- setdiff(lake_pool, exclude_lakes)
  if (!is.null(lakes)) {
    lake_pool <- intersect(lake_pool, lakes)
  }
  lake_pool <- sort(lake_pool)

  results <- purrr::map(lake_pool, function(lk) {
    purrr::imap(
      list(
        Producers = diat_df,
        Consumers = chiro_df
      ),
      function(dca_df, group_label) {
        analyse_lake_gam_timing(
          dca_df = dca_df,
          lake_name = lk,
          group_label = group_label,
          min_n = min_n,
          onset_min_ce = onset_min_ce
        )
      }
    )
  })
  names(results) <- lake_pool

  summary_tbl <- purrr::map_dfr(results, function(x) {
    dplyr::bind_rows(purrr::map(x, "summary"))
  }) %>%
    dplyr::mutate(
      group = factor(group, levels = c("Producers", "Consumers"))
    )

  periods <- purrr::map_dfr(results, function(x) {
    dplyr::bind_rows(purrr::map(x, "periods"))
  }) %>%
    dplyr::mutate(group = factor(group, levels = c("Producers", "Consumers")))

  smooths <- purrr::map_dfr(results, function(x) {
    dplyr::bind_rows(purrr::map(x, "smooth"))
  }) %>%
    dplyr::mutate(group = factor(group, levels = c("Producers", "Consumers")))

  derivs <- purrr::map_dfr(results, function(x) {
    dplyr::bind_rows(purrr::map(x, "deriv"))
  }) %>%
    dplyr::mutate(group = factor(group, levels = c("Producers", "Consumers")))

  comparison <- compare_producer_consumer_timing(summary_tbl, onset_min_ce = onset_min_ce)

  list(
    summary = summary_tbl,
    comparison = comparison,
    periods = periods,
    smooths = smooths,
    derivs = derivs
  )
}

#' Compare producer vs consumer timing within each lake
#'
#' @param timing_summary Output \code{summary} from \code{compute_lake_level_gam_timing}.
#' @param onset_min_ce Age threshold used for onset comparison.
#' @return Tibble with per-lake timing differences and leader labels.
#' @export
compare_producer_consumer_timing <- function(
    timing_summary,
    onset_min_ce = 800) {
  wide <- timing_summary %>%
    dplyr::select(
      lake, group, first_onset_after_min_ce, total_dca_change
    ) %>%
    tidyr::pivot_wider(
      names_from = group,
      values_from = c(first_onset_after_min_ce, total_dca_change),
      names_glue = "{.value}_{group}"
    )

  classify_leader <- function(prod_val, cons_val, tol = 25) {
    if (is.na(prod_val) && is.na(cons_val)) {
      return("Neither detected")
    }
    if (is.na(prod_val)) {
      return("Consumers only")
    }
    if (is.na(cons_val)) {
      return("Producers only")
    }
    diff_val <- prod_val - cons_val
    if (abs(diff_val) <= tol) {
      return("Similar timing")
    }
    if (diff_val < 0) {
      return("Producers earlier")
    }
    "Consumers earlier"
  }

  wide %>%
    dplyr::mutate(
      onset_min_ce = onset_min_ce,
      onset_diff_ce = first_onset_after_min_ce_Producers - first_onset_after_min_ce_Consumers,
      magnitude_diff = total_dca_change_Producers - total_dca_change_Consumers,
      onset_leader = mapply(
        classify_leader,
        first_onset_after_min_ce_Producers,
        first_onset_after_min_ce_Consumers
      )
    ) %>%
    dplyr::arrange(lake)
}

#' Fit a regional HGAM for standardized guild richness
#'
#' @param guild_df Data for one functional guild across lakes.
#' @param min_n Minimum observations required.
#' @param k_global Basis dimension for the global age smooth.
#' @param k_by Basis dimension for lake-specific age smooths.
#' @return Fitted \code{gam} object, or \code{NULL} if fitting fails.
#' @export
fit_regional_guild_richness_gam <- function(
    guild_df,
    min_n = 30L,
    k_global = 30L,
    k_by = 8L) {
  df <- guild_df %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      lake = as.factor(lake),
      age_ce = as.numeric(age_ce),
      standardized_species_richness = as.numeric(standardized_species_richness)
    ) %>%
    dplyr::filter(
      is.finite(age_ce),
      is.finite(standardized_species_richness)
    )

  if (nrow(df) < min_n) {
    return(NULL)
  }

  tryCatch(
    mgcv::gam(
      standardized_species_richness ~
        s(age_ce, bs = "tp", k = k_global) +
        s(age_ce, by = lake, k = k_by, m = 1, bs = "tp") +
        s(lake, bs = "re", k = 8),
      data = df,
      method = "REML"
    ),
    error = function(e) NULL
  )
}

#' Derivative timing for one regional guild-richness HGAM
#'
#' @param guild_df Data for one functional guild across lakes.
#' @param guild_id Guild identifier (e.g. \code{"high_profile"}).
#' @param guild_label Display label for plots/tables.
#' @param trophic_group \code{"Producers"} or \code{"Consumers"}.
#' @param onset_min_ce Ignore significant onsets before this age.
#' @param min_n Minimum observations required to fit the GAM.
#' @return List with \code{summary}, \code{periods}, \code{smooth}, and \code{deriv}.
#' @export
analyse_regional_guild_richness_timing <- function(
    guild_df,
    guild_id,
    guild_label,
    trophic_group,
    onset_min_ce = 800,
    min_n = 30L) {
  n_obs <- guild_df %>%
    dplyr::ungroup() %>%
    dplyr::filter(is.finite(age_ce), is.finite(standardized_species_richness)) %>%
    nrow()

  mod <- fit_regional_guild_richness_gam(guild_df, min_n = min_n)
  if (is.null(mod)) {
    empty <- tibble::tibble(
      guild = guild_id,
      guild_label = guild_label,
      trophic_group = trophic_group,
      n = n_obs,
      first_onset_ce = NA_real_,
      first_onset_after_min_ce = NA_real_,
      first_end_after_min_ce = NA_real_,
      first_direction = NA_character_,
      first_direction_after_min = NA_character_,
      first_peak_ce = NA_real_,
      total_richness_change = NA_real_,
      n_significant_periods = 0L
    )
    return(list(
      summary = empty,
      periods = tibble::tibble(),
      smooth = tibble::tibble(),
      deriv = tibble::tibble()
    ))
  }

  smooth <- gratia::smooth_estimates(mod) %>%
    dplyr::filter(.smooth == "s(age_ce)") %>%
    gratia::add_confint() %>%
    dplyr::transmute(
      guild = guild_id,
      guild_label = guild_label,
      trophic_group = trophic_group,
      age_ce = as.numeric(age_ce),
      estimate = as.numeric(.estimate),
      lower_ci = as.numeric(.lower_ci),
      upper_ci = as.numeric(.upper_ci)
    )

  deriv <- compute_regional_gam_derivatives(mod, guild_label) %>%
    dplyr::mutate(
      guild = guild_id,
      guild_label = guild_label,
      trophic_group = trophic_group
    )

  change <- summarise_gam_change_periods(deriv)
  total_change <- diff(range(smooth$estimate, na.rm = TRUE))

  post <- change$periods %>%
    dplyr::filter(start_ce >= onset_min_ce)

  summary_tbl <- change$summary %>%
    dplyr::transmute(
      guild = guild_id,
      guild_label = guild_label,
      trophic_group = trophic_group,
      n = n_obs,
      first_onset_ce = first_onset_ce,
      first_onset_after_min_ce = if (nrow(post) == 0) NA_real_ else min(post$start_ce),
      first_end_after_min_ce = if (nrow(post) == 0) {
        NA_real_
      } else {
        post$end_ce[which.min(post$start_ce)]
      },
      first_direction = first_direction,
      first_direction_after_min = if (nrow(post) == 0) {
        NA_character_
      } else {
        post$direction[which.min(post$start_ce)]
      },
      first_peak_ce = first_peak_ce,
      total_richness_change = total_change,
      n_significant_periods = n_significant_periods
    )

  periods <- change$periods %>%
    dplyr::mutate(
      guild = guild_id,
      guild_label = guild_label,
      trophic_group = trophic_group
    )

  list(
    summary = summary_tbl,
    periods = periods,
    smooth = smooth,
    deriv = deriv
  )
}

#' Regional guild-richness change timing for all producer/consumer guilds
#'
#' @param prod_df Producer guild richness data.
#' @param cons_df Consumer guild richness data.
#' @param prod_guilds Producer guild ids to include.
#' @param cons_guilds Consumer guild ids to include.
#' @param guild_labels Named vector of display labels.
#' @param exclude_lakes Lakes to exclude.
#' @param onset_min_ce Minimum age for onset reporting.
#' @param min_n Minimum observations per guild GAM.
#' @return List with \code{summary}, \code{periods}, \code{smooths}, and \code{derivs}.
#' @export
compute_regional_guild_richness_timing <- function(
    prod_df,
    cons_df,
    prod_guilds = NULL,
    cons_guilds = NULL,
    guild_labels = NULL,
    exclude_lakes = c("Fogo", "Furnas"),
    onset_min_ce = 800,
    min_n = 30L) {
  prep_df <- function(df) {
    df %>%
      dplyr::ungroup() %>%
      dplyr::filter(
        !.data$lake %in% exclude_lakes,
        is.finite(age_ce),
        is.finite(standardized_species_richness)
      )
  }

  prod_df <- prep_df(prod_df)
  cons_df <- prep_df(cons_df)

  if (is.null(prod_guilds)) {
    prod_guilds <- sort(unique(as.character(prod_df$fgroup)))
  }
  if (is.null(cons_guilds)) {
    cons_guilds <- sort(unique(as.character(cons_df$fgroup)))
  }

  analyse_one <- function(df, guild_id, trophic_group) {
    glab <- if (
      !is.null(guild_labels) && guild_id %in% names(guild_labels)
    ) {
      guild_labels[[guild_id]]
    } else {
      guild_id
    }
    guild_df <- df %>% dplyr::filter(.data$fgroup == guild_id)
    analyse_regional_guild_richness_timing(
      guild_df = guild_df,
      guild_id = guild_id,
      guild_label = glab,
      trophic_group = trophic_group,
      onset_min_ce = onset_min_ce,
      min_n = min_n
    )
  }

  results <- c(
    purrr::map(prod_guilds, ~ analyse_one(prod_df, .x, "Producers")),
    purrr::map(cons_guilds, ~ analyse_one(cons_df, .x, "Consumers"))
  )

  summary_tbl <- purrr::map_dfr(results, "summary") %>%
    dplyr::mutate(
      trophic_group = factor(trophic_group, levels = c("Producers", "Consumers"))
    )

  periods <- purrr::map_dfr(results, "periods") %>%
    dplyr::mutate(
      trophic_group = factor(trophic_group, levels = c("Producers", "Consumers"))
    )

  smooths <- purrr::map_dfr(results, "smooth") %>%
    dplyr::mutate(
      trophic_group = factor(trophic_group, levels = c("Producers", "Consumers"))
    )

  derivs <- purrr::map_dfr(results, "deriv") %>%
    dplyr::mutate(
      trophic_group = factor(trophic_group, levels = c("Producers", "Consumers"))
    )

  list(
    summary = summary_tbl,
    periods = periods,
    smooths = smooths,
    derivs = derivs
  )
}
