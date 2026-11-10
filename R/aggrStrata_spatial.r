# ==============================================================================
# File:      aggrStrata_spatial.r
# Functions: aggrStrataSpatial() [exported]
#
# Purpose:
#   Builds the "strata" dataframe of a domain directly from the sampling frame,
#   given the stratum label of each unit, computing the variance of the target
#   variables in each stratum including the spatial autocorrelation component
#   (parameters fitting, range, kappa). Used by the "spatial" method
#   (optimizeStrataSpatial()).
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
#' @param fitting Fitting of the model(s).
#'   Default is 1.
#' @param range Maximum range for spatial autocorrelation
#' @param kappa Factor used in evaluating spatial autocorrelation.
#'   Default is 3.
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
                              dominio) {

  # 'fitting' is the R^2 of the kriging trend model and appears below as
  # z_z/fitting in the within-stratum variance: fitting = 0 is a division by
  # zero (and any fitting close to 0 correctly forces a very large required
  # sample size, since a model that explains nothing needs near-census
  # sampling to meet any precision target). Reject it here instead of letting
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
  stdev <- function(zz, dist, var, STRATO, dataset, fitting, range, kappa) {
    ind <- which(dataset$STRATO == STRATO)
    # differences
    z_z<-zz[ind,ind]
    # distances
    dist <- dist[ind,ind]
    # variances
    var <- var[ind]

    if (length(ind) > 1) {
      somma_coppie_var <- as.matrix(outer(var,var,"+"))
      spatial_correlation <- (1 - (exp(-kappa*dist/range)))
      # dist = 0 on the diagonal (self-pairs) with range = 0 gives a literal
      # 0/0 = NaN; the correct value there is 0 (zero semivariance with
      # oneself), regardless of range.
      spatial_correlation[dist == 0] <- 0
    }
    if (length(ind) <= 1) {
      somma_coppie_var <- 0
      spatial_correlation <- 0
    }
    # variance in the stratum
    D2 <- z_z/fitting + somma_coppie_var * spatial_correlation
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
  if (!length(fitting) %in% c(1, nvarY) || !length(range) %in% c(1, nvarY))
    stop("aggrStrataSpatial: fitting and range must have length 1 or ", nvarY,
         " (number of Y variables)", call. = FALSE)
  fitting <- rep_len(fitting, nvarY)
  range <- rep_len(range, nvarY)
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
      sd <- stdev(zz,dist,var,j,dataset,fit,rng,kappa)
      stmt <- paste("strato$S",i," <- sd",sep="")
      eval(parse(text=stmt))
    }
    stratatot <- rbind(stratatot, strato)
    }
  return(stratatot)
}