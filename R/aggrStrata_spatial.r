# ==============================================================================
# File:      aggrStrata_spatial.r
# Functions: aggrStrataSpatial() [exported]
#
# Purpose:
#   Builds the "strata" dataframe of a domain directly from the sampling frame,
#   given the stratum label of each unit, computing the variance of the target
#   variables in each stratum including the spatial autocorrelation component
#   (parameters fitting, range, kappa, nugget). Used by KmeansSolutionSpatial().
#   The formula is the same as in buildStrataDFSpatial() (C++), used by the
#   genetic algorithm of optimizeStrataSpatial():
#     D2_ij = (y_i - y_j)^2 / fitting + s2_i + s2_j
#             - 2 (1 - nugget) s_i s_j exp(-kappa d_ij / range),   i != j
#     S2_h  = sum_{i != j} D2_ij / (2 n^2)
#   (until version 1.6 this R version used the covariance (s2_i + s2_j)/2 *
#   exp(...) instead of s_i s_j exp(...), so that K-means and the genetic
#   algorithm evaluated the strata with different formulas).
# ==============================================================================

# ----------------------------------------------------
# Function to aggregate initial atomic strata
# on the basis of the current solution 
# expressed by a vector of integer values
# Version valid for the spatial case
# Author: Giulio Barcaroli
# Date: 6 April 2019
# ----------------------------------------------------
#' @title Builds the "strata" dataframe containing information on target variables Y's distributions in the different strata, starting from a frame where units are spatially correlated.
#'
#' @description
#' This function builds the dataframe "strata" considering as input a given domain in the sampling frame. The variance in each stratum is calculated by considering also the component of spatial autocorrelation.
#'
#' @param dataset This is the name of the dataframe containing the sample data, or the frame data.
#'   It is strictly required that auxiliary information is organised in variables named
#'   as X1, X2, ... , Xm (there should be at least one of them) and the target variables
#'   are denoted by Y1, Y2, ... , Yn.
#'   In addition, in case of sample data, a variable named 'WEIGHT' must be present in the dataframe,
#'   containing   the weigths associated to each sampling unit
#' @param fitting Default is 1, to be used when the frame contains the prediction-error
#'   variances (VAR columns). A value < 1 (R2 of the prediction model) is meant only for
#'   the rescaling correction, with all variances = 0: with positive variances the function
#'   stops (the error variance would be counted twice). If only R2 is known, use fitting = 1
#'   with a constant variance var(predictions) * (1 - R2) / R2. See buildStrataDFSpatial().
#' @param range Range of the exponential correlation function of the prediction errors
#'   (in the units of the coordinates). 0 means no spatial correlation.
#' @param kappa Factor used in evaluating spatial autocorrelation: 1 if range is the
#'   range parameter of an exponential variogram, 3 if it is the practical range.
#'   Any finite value > 0 is accepted; 0 or negative values give an error. Default is 3.
#' @param nugget Relative nugget c0/(c0+c1) of the variogram of the prediction errors,
#'   in [0, 1]. Default is 0.
#' @param vett vector of values indicating how the units in the dataset must
#'   be aggregated in strata.
#' @param dominio Value indicating the domain in the dataset to be processed.
#'
#' @return
#' A dataframe containing strata
#'
#' @author
#' Giulio Barcaroli
#'
#' @examples
#' # synthetic spatial frame: one domain, one target variable with its
#' # prediction variance and the coordinates of each unit
#' set.seed(1234)
#' frame <- data.frame(id = 1:100, X1 = runif(100), Y1 = runif(100, 10, 20),
#'                     var1 = 0.5, lon = runif(100), lat = runif(100),
#'                     domainvalue = 1)
#' # three strata defined on X1
#' strata <- aggrStrataSpatial(dataset = frame,
#'                             fitting = 1,
#'                             range = 0.5,
#'                             kappa = 3,
#'                             nugget = 0.3,
#'                             vett = var.bin(frame$X1, bins = 3),
#'                             dominio = 1)
#' strata
#'
#' @keywords survey
#' @export
aggrStrataSpatial <- function(dataset,
                              fitting=c(1),
                              range=c(0),
                              kappa=3,
                              vett,
                              dominio,
                              nugget=c(0)) {

  # 'fitting' appears below as z_z/fitting in the within-stratum variance:
  # fitting = 0 is a division by zero. Reject it here instead of letting
  # NaN/Inf propagate silently into bethel() and the genetic algorithm.
  if (any(fitting == 0)) {
    stop(
      "aggrStrataSpatial: 'fitting' cannot be 0 (value(s): ",
      paste(fitting, collapse = ", "), ").\n",
      "If your intent is to disable the spatial-correlation adjustment, use ",
      "fitting = 1 together with range = 0 instead of fitting = 0.",
      call. = FALSE
    )
  }

  #---------------------------------------------------------
  # standard deviation calculated with distances
  stdev <- function(zz, dist, var, STRATO, dataset, fitting, range, kappa, nugget) {
    ind <- which(dataset$STRATO == STRATO)
    # differences
    z_z<-zz[ind,ind]
    # distances
    dist <- dist[ind,ind]
    # variances
    var <- var[ind]

    if (length(ind) <= 1) return(0)
    somma_coppie_var <- outer(var, var, "+")
    # covariance of the prediction errors of two distinct units; no spatial
    # correlation with range <= 0 or nugget >= 1
    if (range > 0 && nugget < 1) {
      covar <- (1 - nugget) * sqrt(outer(var, var, "*")) * exp(-kappa * dist / range)
    } else {
      covar <- matrix(0, length(ind), length(ind))
    }
    # model-expected squared differences; self-pairs (i = j) contribute 0
    D2 <- z_z/fitting + somma_coppie_var - 2 * covar
    diag(D2) <- 0
    # variance in the stratum (population variance, divisor n)
    var_strato <- sum(D2) / (2*length(ind)^2)
    # standard deviation
    if (var_strato < 0) var_strato <- 0
    sd_strato <- sqrt(var_strato)
    return(sd_strato)
  }
  #--------------------------------------------------------- 
  colnames(dataset) <- toupper(colnames(dataset))
  dataset <- dataset[dataset$DOMAINVALUE == dominio,]
  dist <- sqrt((outer(dataset$LON,dataset$LON,"-"))^2+(outer(dataset$LAT,dataset$LAT,"-"))^2)
  nvarY <- length(grep("Y", names(dataset)))
  # one value of fitting and range per Y variable: a single value is
  # recycled (otherwise fitting[i], range[i] are NA for i > 1)
  if (!length(fitting) %in% c(1, nvarY) || !length(range) %in% c(1, nvarY) ||
      !length(nugget) %in% c(1, nvarY))
    stop("aggrStrataSpatial: fitting, range and nugget must have length 1 or ", nvarY,
         " (number of Y variables)", call. = FALSE)
  if (length(kappa) != 1 || !is.finite(kappa) || kappa <= 0)
    stop("'kappa' must be a single finite value > 0", call. = FALSE)
  if (any(is.na(nugget)) || any(nugget < 0 | nugget > 1))
    stop("aggrStrataSpatial: 'nugget' must be in [0, 1] (relative nugget c0/(c0+c1))",
         call. = FALSE)
  fitting <- rep_len(fitting, nvarY)
  range <- rep_len(range, nvarY)
  nugget <- rep_len(nugget, nvarY)
  for (i in 1:nvarY) {
    vi <- dataset[[paste0("VAR", i)]]
    if (fitting[i] < 1 && !is.null(vi) && any(vi > 0, na.rm = TRUE))
      stop("fitting < 1 cannot be used together with positive prediction-error variances: ",
           "the error variance would be counted twice. Either omit fitting (default 1) ",
           "and keep the prediction-error variances in the frame, or, if only the R2 of ",
           "the prediction model is known, use fitting = 1 with a constant variance ",
           "equal to var(predictions) * (1 - R2) / R2.", call. = FALSE)
  }
  stratatot <- NULL
  dataset$STRATO <- as.factor(vett)
  for (j in as.numeric(levels(dataset$STRATO))) {
    strato <- NULL
    strato$stratum <- j
    strato$N <- nrow(dataset[dataset$STRATO == j,])
    strato$COST <- 1
    strato$CENS <- 0
    strato$DOM1 <- dominio
    strato <- as.data.frame(strato, stringsAsFactors = TRUE)
    rng <- NULL
    fit <- NULL
    zz <- NULL
    for (i in 1:nvarY) {
      stmt <- paste("var <- dataset$VAR",i,sep="")
      eval(parse(text = stmt))
      stmt <- paste("strato$M", i, " <- mean(dataset$Y",i,"[dataset$STRATO == ",j,"])",sep = "")
      eval(parse(text = stmt))
      stmt <- paste("zz <- outer(dataset$Y",i,",dataset$Y",i,",'-')^2",sep="")
      eval(parse(text = stmt))
      stmt <- paste("rng <- range[",i,"]",sep="")
      eval(parse(text = stmt))
      stmt <- paste("fit <- fitting[",i,"]",sep="")
      eval(parse(text = stmt))
      sd <- stdev(zz,dist,var,j,dataset,fit,rng,kappa,nugget[i])
      stmt <- paste("strato$S",i," <- sd",sep="")
      eval(parse(text=stmt))
    }
    stratatot <- rbind(stratatot, strato)
    }
  return(stratatot)
}