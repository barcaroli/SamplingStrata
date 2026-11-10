# ==============================================================================
# File:    SamplingStrata-package.R
# Purpose: package-level roxygen directives used to generate the NAMESPACE
#          file: registration of the compiled code (Rcpp) and the functions
#          imported from other packages. This file contains no executable
#          code and generates no documentation (@noRd).
# ==============================================================================

#' @useDynLib SamplingStrata, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom grDevices dev.off pdf png rainbow
#' @importFrom graphics abline barplot boxplot close.screen legend lines par
#'   plot points screen split.screen title
#' @importFrom stats aggregate cov dnorm kmeans runif sd var
#' @importFrom utils globalVariables install.packages read.delim
#'   setTxtProgressBar txtProgressBar write.table
#' @importFrom memoise memoise
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel clusterExport detectCores makePSOCKcluster stopCluster
#' @importFrom pbapply pblapply
#' @importFrom SamplingBigData lpm2_kdtree
#' @importFrom glue glue
#' @importFrom methods is
#' @noRd
NULL
