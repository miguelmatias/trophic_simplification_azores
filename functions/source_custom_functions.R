
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


# INPUTS

#' Get optimum AMD clusters
#'
#' @param
#' @returns
#' @seealso
#' @examples
#' @export
optAMDclusters <- function(.data, .iterations = 100, .min_clusts = 2, .num_groups = 12) {

  .maxpms <- c() # set vector to collect permutations
  
  for (k in 1:.iterations) { # loop over set number of iterations
    
    .maxpm <- c() # set partial permutation vector
    
    for (i in .min_clusts:.num_groups) { # loop range of min to max clusters
      
      .clusts_k <- cmeans(.data, i, 20, verbose = F, method = "cmeans", m = 2) # fuzzy clustering per iteration k
      
      .prob <- .clusts_k$membership # membership probability of each group for each sample
      
      # Selecciona la probabilidad de pertenencia de cada muestra al grupo al que fue asignada en "cluster".
      .maxprob <- apply(.prob, 1, max)
      
      .mpm <- mean(.maxprob) - 1 / i # Calcula la media de las probabilidades max. calculadas en el paso anterior.
      .maxpm <- cbind(.maxpm, .mpm) # crea una BD con las Mpm de cada uno de los an?lisis, con distinto n? de grupos, de .min_clusts a .num_groups.
      n <- ncol(.maxpm)
      colnames(.maxpm)[n] <- paste("C", i, sep = "") # nombra el n? de clusters al que corresponde cada valor de CRc
    }
    
    .maxpms <- rbind(.maxpms, .maxpm) # crea un BD con los valores de .maxpm de cada itr
    
    #if (i == .num_groups) print(k) # para que vaya indicando, mientras corre el script, por qu? n? de itr va
  }
  
  # set function to pick maximum
  colMax <- function(.col_data) {
    apply(.col_data, MARGIN = c(2), max)
  }
  
  # results
  .maxpmean <- colMax(.maxpms) # calcula el max. de todos los valores de CRs obtenidos para cada cada n? de clusters
  .results <- bind_cols(.clusters = c(.min_clusts:.num_groups), .maxpmean = .maxpmean)
  
  # get max AMD value
  .maxAMD <- round(max(.results$.maxpmean), 2)
  
  # get number of clusters
  .opt_clust_numb <- .results$.clusters[which(.results$.maxpmean == max(.results$.maxpmean))]
  
  # collect data as results list
  .output <- list(optimal = .opt_clust_numb)
  
  return(.output$optimal)
}

optAMDclusters_2 <- function(.data, .iterations = 100, .min_clusts = 2, .num_groups = 12) {
  
  # Helper function: Calculates the maximum mean probability adjusted by the number of clusters
  calc_mpm <- function(prob_matrix, num_clusters) {
    .maxprob <- apply(prob_matrix, 1, max) # Max probability per sample
    mean(.maxprob) - 1 / num_clusters      # Adjusted mean
  }
  
  # Preallocate the matrix to store results for each iteration
  .maxpms <- matrix(0, nrow = .iterations, ncol = (.num_groups - .min_clusts + 1))
  colnames(.maxpms) <- paste0("C", .min_clusts:.num_groups)
  
  # Main iteration loop
  for (k in seq_len(.iterations)) {
    # Compute the mpm for each number of clusters in the range
    .maxpm <- sapply(.min_clusts:.num_groups, function(i) {
      .clusts_k <- cmeans(.data, i, 20, verbose = FALSE, method = "cmeans", m = 2)
      calc_mpm(.clusts_k$membership, i)
    })
    .maxpms[k, ] <- .maxpm
  }
  
  # Compute column maxima
  .maxpmean <- apply(.maxpms, 2, max)
  
  # Identify optimal number of clusters
  .results <- tibble::tibble(
    .clusters = .min_clusts:.num_groups,
    .maxpmean = .maxpmean
  )
  .maxAMD <- round(max(.results$.maxpmean), 2)
  .opt_clust_numb <- .results$.clusters[which.max(.results$.maxpmean)]
  
  # Return optimal number of clusters
  return(.opt_clust_numb)
}

# INPUTS

#' Plot AMD curve with optimum AMD clusters
#'
#' @param
#' @returns
#' @seealso
#' @examples
#' @export
PlotoptAMDclusters <- function(.data, .iterations = 100, .min_clusts = 2, .num_groups = 12) {

  .maxpms <- c() # set vector to collect permutations
  
  for (k in 1:.iterations) { # loop over set number of iterations
    
    .maxpm <- c() # set partial permutation vector
    
    for (i in .min_clusts:.num_groups) { # loop range of min to max clusters
      
      .clusts_k <- cmeans(.data, i, 20, verbose = F, method = "cmeans", m = 2) # fuzzy clustering per iteration k
      
      .prob <- .clusts_k$membership # membership probability of each group for each sample
      
      # Selecciona la probabilidad de pertenencia de cada muestra al grupo al que fue asignada en "cluster".
      .maxprob <- apply(.prob, 1, max)
      
      .mpm <- mean(.maxprob) - 1 / i # Calcula la media de las probabilidades max. calculadas en el paso anterior.
      .maxpm <- cbind(.maxpm, .mpm) # crea una BD con las Mpm de cada uno de los an?lisis, con distinto n? de grupos, de .min_clusts a .num_groups.
      n <- ncol(.maxpm)
      colnames(.maxpm)[n] <- paste("C", i, sep = "") # nombra el n? de clusters al que corresponde cada valor de CRc
    }
    
    .maxpms <- rbind(.maxpms, .maxpm) # crea un BD con los valores de .maxpm de cada itr
    
    #if (i == .num_groups) print(k) # para que vaya indicando, mientras corre el script, por qu? n? de itr va
  }
  
  # set function to pick maximum
  colMax <- function(.col_data) {
    apply(.col_data, MARGIN = c(2), max)
  }
  
  # results
  .maxpmean <- colMax(.maxpms) # calcula el max. de todos los valores de CRs obtenidos para cada cada n? de clusters
  .results <- bind_cols(.clusters = c(.min_clusts:.num_groups), .maxpmean = .maxpmean)
  
  # get max AMD value
  .maxAMD <- round(max(.results$.maxpmean), 2)
  
  # get number of clusters
  .opt_clust_numb <- .results$.clusters[which(.results$.maxpmean == max(.results$.maxpmean))]
  
  # collect data as results list
  .output <- list(
    optimal = .opt_clust_numb,
    results = .results
    )
  
  return(.output)
}


util_make_trend <-
  function(data_source,
           sel_method = c("linear", "non_linear")) {
    util_check_class("data_source", "data.frame")
    
    util_check_col_names("data_source", c("ROC", "Age"))
    
    util_check_class("sel_method", "character")
    
    util_check_vector_values("sel_method", c("linear", "non_linear"))
    
    sel_method <- match.arg(sel_method)
    
    if (
      sel_method == "non_linear"
    ) {
      res <-
        mgcv::predict.gam(
          mgcv::gam(
            ROC ~ s(Age, k = 3),
            data = data_source,
            family = mgcv::Tweedie(p = 2),
            method = "REML"
          ),
          type = "response"
        )
    } else {
      res <-
        stats::predict.glm(
          stats::glm(ROC ~ Age,
                     data = data_source,
                     family = mgcv::Tweedie(p = 2)
          ),
          type = "response"
        )
    }
    
    return(res)
  }

#Rate of change functions

util_check_class <-
  function(data_source, sel_class) {
    parent_frame <- sys.parent()
    
    parent_env <- sys.frame(which = parent_frame)
    
    data_source_class <- class(get(data_source, envir = parent_env))
    
    assertthat::assert_that(
      any(data_source_class %in% sel_class),
      msg = paste0(
        "'", data_source, "' must be one of the following: ",
        util_paste_as_vector(sel_class)
      )
    )
  }

util_paste_as_vector <-
  function(var_list, sep = "'") {
    paste(
      paste0(sep, var_list, sep),
      collapse = ", "
    ) %>%
      return()
  }

util_check_col_names <-
  function(data_source, var_list) {
    parent_frame <- sys.parent()
    
    parent_env <- sys.frame(which = parent_frame)
    
    data_source_obj <- get(data_source, envir = parent_env)
    
    util_check_class("data_source_obj", "data.frame")
    
    util_check_class("var_list", "character")
    
    assertthat::assert_that(
      all(var_list %in% names(data_source_obj)),
      msg = paste0(
        "'", data_source, "' must contains following columns: ",
        util_paste_as_vector(var_list)
      )
    )
  }

util_check_vector_values <-
  function(data_source, var_list) {
    parent_frame <- sys.parent()
    
    parent_env <- sys.frame(which = parent_frame)
    
    data_source_obj <- get(data_source, envir = parent_env)
    
    util_check_class("data_source_obj", "character")
    
    util_check_class("var_list", "character")
    
    assertthat::assert_that(
      any(var_list %in% data_source_obj),
      msg = paste0(
        "'", data_source, "' must contains one of the following values: ",
        util_paste_as_vector(var_list)
      )
    )
  }

# CUSTOM GLASSO RELATED FUNCTIONS

#########################################
# Function for sampleSize:
sampleSize_pairwise <- function(data, type = c( "pairwise_average","maximum","minimum","pairwise_maximum",
                                                "pairwise_minimum")){
  type <- match.arg(type)
  
  if (type == "maximum"){
    
    sampleSize <- sum(apply(data,1,function(x)!all(is.na(x))))
    
  } else if (type == "minimum"){
    
    sampleSize <- sum(apply(data,1,function(x)!any(is.na(x))))
    
  } else {
    # Matrix with NAs:
    xmat <- as.matrix(!is.na(data))
    # product:
    misMatrix <- t(xmat) %*% xmat
    
    if (type == "pairwise_maximum"){
      sampleSize <- max(misMatrix)
    } else  if (type == "pairwise_minimum"){
      sampleSize <- min(misMatrix)
    } else  if (type == "pairwise_average"){
      sampleSize <- mean(misMatrix)
    }
  }
  
  return(sampleSize)
}

### EBIC GLASSO ESTIMATOR ###
bootnet_EBICglasso <- function(
    data, # Dataset used
    tuning = 0.5, # tuning parameter
    corMethod = c("cor","cov","cor_auto","npn","spearman"), # Correlation method
    missing = c("pairwise","listwise","fiml","stop"),
    sampleSize = c("pairwise_average","maximum","minimum","pairwise_maximum",
                   "pairwise_minimum"), # Sample size when using missing = "pairwise"
    verbose = TRUE,
    corArgs = list(), # Extra arguments to the correlation function
    refit = FALSE,
    principalDirection = FALSE,
    lambda.min.ratio = 0.01,
    nlambda = 100,
    threshold = FALSE,
    unlock = FALSE,
    nonPositiveDefinite = c("stop","continue"),
    transform = c("none","rank","quantile"),
    ...){
  
  transform <- match.arg(transform)
  if (transform == "rank"){
    data <- rank_transformation(data)
  } else if (transform == "quantile"){
    data <- quantile_transformation(data)
  }
  
  nonPositiveDefinite <- match.arg(nonPositiveDefinite)
  if (!unlock){
    stop("You are using an internal estimator function without using 'estimateNetwork'. This function is only intended to be used from within 'estimateNetwork' and will not run now. To force manual use of this function (not recommended), use unlock = TRUE.")  
  }
  
  # Check arguments:
  corMethod <- match.arg(corMethod)
  missing <- match.arg(missing)
  # sampleSize <- match.arg(sampleSize)
  
  # Message:
  if (verbose){
    msg <- "Estimating Network. Using package::function:"  
    msg <- paste0(msg,"\n  - qgraph::EBICglasso for EBIC model selection\n    - using glasso::glasso")
    if (corMethod == "cor_auto"){
      msg <- paste0(msg,"\n  - qgraph::cor_auto for correlation computation\n    - using lavaan::lavCor")
    }
    if (corMethod == "npn"){
      msg <- paste0(msg,"\n  - huge::huge.npn for nonparanormal transformation")
    }
    # msg <- paste0(msg,"\n\nPlease reference accordingly\n")
    message(msg)
  }
  
  
  # First test if data is a data frame:
  if (!(is.data.frame(data) || is.matrix(data))){
    stop("'data' argument must be a data frame")
  }
  
  # If matrix coerce to data frame:
  if (is.matrix(data)){
    data <- as.data.frame(data)
  }
  
  # Obtain info from data:
  N <- ncol(data)
  Np <- nrow(data)
  
  # Check missing:
  if (missing == "stop"){
    if (any(is.na(data))){
      stop("Missing data detected and missing = 'stop'")
    }
  }
  
  # Correlate data:
  corMat <- bootnet_correlate2(data = data, corMethod =  corMethod, 
                               corArgs = corArgs, missing = missing,
                               verbose = verbose,nonPositiveDefinite=nonPositiveDefinite)
  
  
  
  
  # Sample size:
  if (missing == "listwise"){
    sampleSize <- nrow(na.omit(data))
  } else{
    sampleSize <- sampleSize_pairwise(data, sampleSize)
  } 
  
  # Principal direction:
  if (principalDirection){
    corMat <- principalDirection(corMat)
  }
  
  # Estimate network:
  Results <- qgraph::EBICglasso(corMat,
                                n =  sampleSize, 
                                gamma = tuning,
                                returnAllResults = TRUE,
                                refit = refit,
                                lambda.min.ratio=lambda.min.ratio,
                                nlambda = nlambda,
                                threshold=threshold,
                                ...)
  
  # Return:
  return(list(graph=Results$optnet,results=Results))
}



# Function that checks input and returns the functions:
checkInput2 <- function(
    default = c("none", "EBICglasso","ggmModSelect", "pcor","IsingFit","IsingSampler", "huge","adalasso","mgm","relimp", 
                "cor","TMFG","ggmModSelect","LoGo","graphicalVAR","piecewiseIsing","SVAR_lavaan",
                "GGMncv"),
    fun, # Estimator function
    # prepFun, # Fun to produce the correlation or covariance matrix
    # prepArgs, # list with arguments for the correlation function
    # estFun, # function that results in a network
    # estArgs, # arguments sent to the graph estimation function (if missing automatically sample size is included)
    # graphFun, # set to identity if missing
    # graphArgs, # Set to null if missing
    # intFun, # Set to null if missing
    # intArgs, # Set to null if missing
    # nSample,
    verbose=TRUE,
    # construct = c("default","function","arguments"),
    .dots = list(),
    ... # Arguments to the estimator function
){
  construct <- "function"
  if (default[[1]]=="glasso") default <- "EBICglasso"
  if (default[[1]]=="IsingSampler") default <- "IsingSampler"
  default <- match.arg(default)
  # construct <- match.arg(construct)
  
  ### DEFAULT OPTIONS ###
  if (missing(fun)){
    fun <- NULL
  }
  
  
  
  # Stop if not compatible:
  dots <- c(.dots,list(...))
  
  # gather names:
  argNames <- character(0)
  # 
  # if (!missing(prepFun)){
  #   argNames <- c(argNames,"prepFun")
  # }
  # if (!missing(prepArgs)){
  #   argNames <- c(argNames,"prepArgs")
  # }
  # if (!missing(estFun)){
  #   argNames <- c(argNames,"estFun")
  # }
  # if (!missing(estArgs)){
  #   argNames <- c(argNames,"estArgs")
  # }
  # if (!missing(graphFun)){
  #   argNames <- c(argNames,"graphFun")
  # }
  # if (!missing(graphArgs)){
  #   argNames <- c(argNames,"graphArgs")
  # }
  # if (!missing(intFun)){
  #   argNames <- c(argNames,"intFun")
  # }
  # if (!missing(intArgs)){
  #   argNames <- c(argNames,"intArgs")
  # }
  # 
  # # Not compatible if construct is used:
  # if (length(dots) > 0 && construct == "arguments"){
  #   
  #   stop(paste0("Ambiguous argument specification. Old functonality is used (construct = 'arguments') in combination with new functionality arguments (implying construct = 'function'): ",
  #               paste0("'",names(dots),"'",collapse="; "),". These arguments are NOT compatible!"))
  #   
  # }
  # 
  # # relimp not compatable with old:
  # if (construct == "arguments" & default == "relimp"){
  #   stop("default = 'relimp' not supported with old bootnet style (construct = 'arguments')")
  #   
  # }
  # 
  # if (length(argNames) > 0 && construct == "function"){
  #   
  #   stop(paste0("Ambiguous argument specification. New functonality is used (construct = 'function') in combination with old functionality arguments (implying construct = 'arguments'): ",
  #               paste0("'",argNames,"'",collapse="; "),". These arguments are NOT compatible!"))
  #   
  # }
  #   
  # # not compatible if both dots are used and arguments are used:
  # if (length(argNames) > 0 & length(dots) > 0){
  # 
  #   stop(paste0("Ambiguous argument specification. Both old functionality arguments are used, compatible with construct = 'arguments': ",
  #               paste0("'",argNames,"'",collapse="; "),", as well as new functionality arguments are used, compatible with construct = 'function': ",
  #               paste0("'",names(dots),"'",collapse="; "),". These two types of arguments are NOT compatible!"))
  #   
  # }
  # 
  # 
  # # Check to construct via function or to construct via arguments:
  # # if no default and no fun, use arguments:
  # if (construct == "default"){
  #   construct <- "function"
  #   
  #   if (default == "none" && is.null(fun)){
  #     construct <- "arguments"
  #   }
  #   
  #   # If fun is missing, default is not none and one argument is not missing, use arguments (backward competability):
  #   if (default != "none" && is.null(fun) && (!missing(prepFun) | !missing(prepArgs) | !missing(estFun) | !missing(estArgs))){
  #     construct <- "arguments"
  #   }
  # }
  # 
  # # Check if arguments are not missing:
  # if (default == "none" && construct == "arguments"){
  #   if (missing(prepFun) | missing(prepArgs) | missing(estFun) | missing(estArgs)){
  #     stop("If 'default' is not set and 'fun' is missing, 'prepFun', 'prepArgs', 'estFun' and 'estArgs' may not be missing.")
  #   }
  # }
  
  ### Construct estimator function via function:
  if (construct == "function"){
    # Arguments:
    Args <- dots
    # 
    # # Warn user that arguments are ignored:
    # if (!missing(prepFun)){
    #   warning("'prepFun' argument is ignored as a function is used as arguments. To use 'prepFun', please set construct = 'arguments'")
    # }
    # if (!missing(prepArgs)){
    #   warning("'prepArgs' argument is ignored as a function is used as arguments. To use 'prepArgs', please set construct = 'arguments'")
    # }
    # if (!missing(estFun)){
    #   warning("'estFun' argument is ignored as a function is used as arguments. To use 'estFun', please set construct = 'arguments'")
    # }
    # if (!missing(estArgs)){
    #   warning("'estArgs' argument is ignored as a function is used as arguments. To use 'estArgs', please set construct = 'arguments'")
    # }
    # if (!missing(graphFun)){
    #   warning("'graphFun' argument is ignored as a function is used as arguments. To use 'graphFun', please set construct = 'arguments'")
    # }
    # if (!missing(graphArgs)){
    #   warning("'graphArgs' argument is ignored as a function is used as arguments. To use 'graphArgs', please set construct = 'arguments'")
    # }
    # if (!missing(intFun)){
    #   warning("'intFun' argument is ignored as a function is used as arguments. To use 'intFun', please set construct = 'arguments'")
    # }
    # if (!missing(intArgs)){
    #   warning("'intArgs' argument is ignored as a function is used as arguments. To use 'intArgs', please set construct = 'arguments'")
    # }
    # 
    # per default:
    if (default == "none"){
      Function <- fun
    } else if (default == "EBICglasso"){
      Function <- bootnet_EBICglasso
    } else if (default == "ggmModSelect"){
      Function <- bootnet_ggmModSelect
    } else if (default == "IsingFit"){
      Function <- bootnet_IsingFit
    } else if (default == "IsingSampler"){
      Function <- bootnet_IsingSampler
    } else if (default == "pcor"){
      Function <- bootnet_pcor
    } else if (default == "cor"){
      Function <- bootnet_cor
    } else if (default == "adalasso"){
      Function <- bootnet_adalasso
    } else if (default == "huge"){
      Function <- bootnet_huge
    } else if (default == "mgm"){
      Function <- bootnet_mgm
    } else if (default == "relimp"){
      Function <- bootnet_relimp
    } else if (default == "TMFG"){
      Function <- bootnet_TMFG
    } else if (default == "LoGo"){
      Function <- bootnet_LoGo
    } else if (default == "graphicalVAR"){
      Function <- bootnet_graphicalVAR  
    } else if (default == "piecewiseIsing"){
      Function <- bootnet_piecewiseIsing
    } else if (default == "SVAR_lavaan"){
      Function <- bootnet_SVAR_lavaan  
    } else if (default == "GGMncv"){
      Function <- bootnet_GGMncv 
    } else stop("Currently not supported.")
    
  }
  
  # Output:
  Output <- list(
    data = data,
    default = default,
    estimator = Function,
    arguments = Args
  )
  
  return(Output)
}

# Function for correlation/covariance:
bootnet_correlate2 <- function(data, corMethod =  c("cor","cor_auto","cov","npn","spearman"), 
                               corArgs = list(), missing = c("pairwise","listwise","fiml","stop"),
                               verbose = TRUE, nonPositiveDefinite = "continue",
                               transform = c("none","rank","quantile")){
  transform <- match.arg(transform)
  if (transform == "rank"){
    data <- rank_transformation(data)
  } else if (transform == "quantile"){
    data <- quantile_transformation(data)
  }
  
  nonPositiveDefinite <- match.arg(nonPositiveDefinite)
  corMethod <- match.arg(corMethod)
  missing <- match.arg(missing)
  
  # Correlate data:
  # npn:
  if (corMethod == "npn"){
    data <- huge::huge.npn(data)
    corMethod <- "cor"
  }
  
  # cor_auto:
  if (corMethod == "cor_auto"){
    args <- list(data=data,missing=missing,verbose=verbose)
    if (length(corArgs) > 0){
      for (i in seq_along(corArgs)){
        args[[names(corArgs)[[i]]]] <- corArgs[[i]]
      }
    }
    corMat <- do.call(qgraph::cor_auto,args)
  } else if (corMethod%in%c("cor","cov","spearman")){
    # Normal correlations
    if (missing == "fiml"){
      stop("missing = 'fiml' only supported with corMethod = 'cor_auto'")
    }
    use <- switch(missing,
                  pairwise = "pairwise.complete.obs",
                  listwise = "complete.obs")
    
    args <- list(x=data,use=use)
    if (length(corArgs) > 0){
      for (i in seq_along(corArgs)){
        args[[names(corArgs)[[i]]]] <- corArgs[[i]]
      }
    }
    if (corMethod == "spearman"){
      args[["method"]] <- "spearman"
      corMethod <- "cor"
    }
    
    corMat <- do.call(corMethod,args)
    corMat <- corpcor::make.positive.definite(corMat)
  } else stop ("Correlation method is not supported.")
  
  if (nonPositiveDefinite == "stop"){
    if (!all(eigen(corMat)$values > 0)){
      stop("Correlation matrix is not positive definite.")
    }    
  }
  
  return(corMat)
}

# This function takes data as input and produced a network. It is used inside bootnet:
estimateNetwork2 <- function(
    data,
    default = c("none", "EBICglasso", "pcor","IsingFit","IsingSampler", "huge","adalasso","mgm","relimp", "cor","TMFG",
                "ggmModSelect", "LoGo","graphicalVAR", "piecewiseIsing","SVAR_lavaan",
                "GGMncv"),
    fun, # A function that takes data and returns a network or list entitled "graph" and "thresholds". optional.
    # prepFun, # Fun to produce the correlation or covariance matrix
    # prepArgs, # list with arguments for the correlation function
    # estFun, # function that results in a network
    # estArgs, # arguments sent to the graph estimation function (if missing automatically sample size is included)
    # graphFun, # set to identity if missing
    # graphArgs, # Set to null if missing
    # intFun, # Set to null if missing
    # intArgs, # Set to null if missing
    labels, # if missing taken from colnames
    verbose = TRUE, # Dummy used in cor_auto and in the future in other functions. Set to FALSE in bootnet
    # construct = c("default","function","arguments"),
    .dots = list(),
    weighted = TRUE,
    signed = TRUE,
    directed,
    datatype,
    checkNumeric = FALSE,
    # plot = TRUE, # Plot the network?
    ..., # Arguments to the 'fun' function
    .input, # Skips most of first steps if supplied
    memorysaver = FALSE # If set to FALSE data, estimator and results are not stored.
){
  construct <- "function"
  # Borsboom easter egg:
  if (default[1] == "Borsboom") return(42)
  
  if (default[[1]]=="glasso") default <- "EBICglasso"
  default <- match.arg(default)
  
  # datatype test:
  if (missing(datatype)){
    if (is(data,"tsData")){
      datatype <- "graphicalVAR"
    } else {
      datatype <- "normal"
    }
    
  }
  
  if (!datatype%in% c("normal","graphicalVAR")){
    stop("Only datatypes 'normal' and 'graphicalVAR' currently supported.")
  }
  #   
  # If NAs and default can't handle, stop:
  # if (any(is.na(data)) && default %in% c("huge","adalasso")){
  #   stop(paste0("Missing data not supported for default set '",default,"'. Try using na.omit(data)."))
  # }
  
  # First test if data is a data frame:
  if (datatype == "normal" && !(is.data.frame(data) || is.matrix(data))){
    stop("'data' argument must be a data frame")
  }
  
  # If matrix coerce to data frame:
  if (datatype == "normal" && is.matrix(data)){
    data <- as.data.frame(data)
  }
  
  if (missing(directed)){
    if (default == "graphicalVAR"){
      directed <- list(contemporaneous = FALSE, temporal = TRUE)
    } else  if (default == "SVAR_lavaan"){
      directed <- list(contemporaneous = TRUE, temporal = TRUE)
    } else if (!default %in% c("relimp","DAG")){
      directed <- FALSE 
    } else {
      directed <- TRUE
    }
  }
  
  if (datatype == "normal"){
    N <- ncol(data)
    Np <- nrow(data)   
    if (missing(labels)){
      labels <- colnames(data)
    }
    
    if (checkNumeric){
      # Check and remove any variable that is not ordered, integer or numeric:
      goodColumns <- sapply(data, function(x) is.numeric(x) | is.ordered(x) | is.integer(x))
      
      if (!all(goodColumns)){
        if (verbose){
          warning(paste0("Removing non-numeric columns: ",paste(which(!goodColumns),collapse="; ")))
        }
        data <- data[,goodColumns,drop=FALSE]
      }
    }
    
    
  } else if (datatype == "graphicalVAR"){
    N <- length(data$vars)
    Np <- nrow(data$data_c)
    if (missing(labels)){
      labels <- data$vars
    }
  }
  
  # Compute estimator:
  if (missing(.input)){
    .input <- checkInput2(
      default = default,
      fun = fun,
      # prepFun = prepFun, # Fun to produce the correlation or covariance matrix
      # prepArgs = prepArgs, # list with arguments for the correlation function
      # estFun=estFun, # function that results in a network
      # estArgs=estArgs, # arguments sent to the graph estimation function (if missing automatically sample size is included)
      # graphFun=graphFun, # set to identity if missing
      # graphArgs=graphArgs, # Set to null if missing
      # intFun=intFun, # Set to null if missing
      # intArgs=intArgs, # Set to null if missing
      # sampleSize = Np,
      # construct=construct,
      verbose=verbose,
      .dots=.dots,
      ...
    )
  }
  
  
  # Add verbose:
  # Every estimator must have argument verbose:
  if ("verbose" %in% names(formals(.input$estimator))){
    .input$arguments$verbose <- verbose
  }
  
  # Unlock function:
  # Every estimator must have argument verbose:
  if ("unlock" %in% names(formals(.input$estimator))){
    .input$arguments$unlock <- TRUE
  }
  
  # Compute network:
  Result <- do.call(.input$estimator, c(list(data),.input$arguments))
  
  if (!is.list(Result)){
    sampleGraph <- Result
    intercepts <- NULL
    output <- Result
    nNode <- ncol(Result)
  } else if (is.list(Result$graph)){
    sampleGraph <- Result$graph
    intercepts <- Result$intercepts
    output <- Result$results
    nNode <- ncol(Result$graph[[1]])
  } else {
    sampleGraph <- Result$graph
    intercepts <- Result$intercepts
    output <- Result$results
    nNode <- ncol(Result$graph)
  }
  
  
  if (!is.matrix(sampleGraph)){
    if (is.list(sampleGraph)){
      if (!is.matrix(sampleGraph[[1]])){
        stop("Estimated result is not a list of matrices encoding networks.")   
      }
    } else {
      stop("Estimated result is not a matrix encoding a network.")
    }
  }
  
  # Special data?
  if (!is.list(Result) || is.null(Result$specialData)){
    outdata <- data
    datatype <- "normal"
  } else {
    outdata <- Result$specialData$data
    datatype <- Result$specialData$type
  }
  
  
  sampleResult <- list(
    graph = sampleGraph,
    intercepts = intercepts,
    results = output,
    labels = labels,
    nNode = nNode,
    nPerson = Np,
    estimator = .input$estimator,
    arguments = .input$arguments,
    data = outdata,
    datatype = datatype,
    default = default,
    weighted = weighted,
    signed = signed,
    directed=directed,
    .input = .input,
    thresholded = FALSE
  )
  class(sampleResult) <- c("bootnetResult", "list")
  
  if (default == "graphicalVAR"){
    sampleResult$labels <- output$data$vars
  }
  if (default == "SVAR_lavaan"){
    sampleResult$labels  <- outdata$vars
  }
  
  # Memory save:
  if(memorysaver)
  {
    sampleResult$results <- NA
    sampleResult$estimator <- NA
    sampleResult$data <- NA
    sampleResult$.input <- NA
  }
  
  # Plot?
  #   if (plot){
  #     plot(sampleResult,labels=labels,layout = "spring", ...)
  #   }
  
  # Return network:
  return(sampleResult)
}

