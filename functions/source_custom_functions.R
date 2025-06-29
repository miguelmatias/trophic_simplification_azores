
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
subsample_time_series <- function(data, x, y, binwidth, n_samples_per_bin) {
  # Ensure reproducibility
  set.seed(123)
  
  # Create a bin column based on the binwidth
  data$bin <- cut(data[[x]], breaks = seq(min(data[[x]]), max(data[[x]]), by = binwidth), include.lowest = TRUE)
  
  # Randomly subsample the same number of points from each bin
  balanced_data <- data %>%
    group_by(bin) %>%
    sample_n(size = min(n(), n_samples_per_bin), replace = FALSE)
  
  # Return the subsampled data
  return(balanced_data)
}

#' Fix slumps / add decimals
#'
#' @param
#' @returns
#' @seealso
#' @examples
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

#' Get brokenstick values for each dataset
#'
#' @param
#' @returns
#' @seealso
#' @examples
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

#' Get lake names from dataframes within lists
#'
#' @param
#' @returns
#' @seealso
#' @examples
#' @export
get_names_list <- function(.lista) {
  .names <- lapply(.lista, function(x) unique(x[, "lake"]))
  return(.names)
}

#' Get AMD clusters
#'
#' @param
#' @returns
#' @seealso
#' @examples
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
