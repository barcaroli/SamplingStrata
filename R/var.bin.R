# ==============================================================================
# File:      var.bin.R
# Functions: var.bin() [exported]
#
# Purpose:
#   Transforms a continuous variable into an ordinal categorical one with a
#   given number of classes, by k-means clustering on the variable, so that it
#   can be used as stratification variable by the "atomic" method.
# ==============================================================================

# ----------------------------------------------------
# Function to transform cotinuous stratification
# variables into ordered factors by using the
# k-means clustered method
# Author: Giulio Barcaroli
# Date: 4 January 2012
# ----------------------------------------------------
#' @title Allows to transform a continuous variable into a categorical ordinal one by applying a modified version of the k-means clustering function in the 'stats' package.
#'
#' @description
#' The optimization of a frame stratification is applicable only in presence of all categorical
#' auxiliary variables in the frame. If one or more continuous auxiliary variables are in the frame,
#' it is necessary to pre-process in order to convert them into categorical (ordinal) variables.
#' The applied method is the "k-means" clustering method contained in the in "stats" package.
#' This function ensures that the final result is in an ordered categorical variable.
#'
#' @param x Continuous variable to be transformed into a categorical one
#' @param bins Number of values of the resulting categorical variable
#' @param iter.max Maximum number of iterations of the clustering algorithm
#'
#' @return
#' Binned variable
#'
#' @examples
#' \donttest{
#' library(SamplingStrata)
#' data(swissmunicipalities)
#' data(swissframe)
#' swissframe$X1 <- var.bin(swissmunicipalities$POPTOT,bins = 18)
#' table(swissframe$X1)
#' tapply(swissmunicipalities$POPTOT,swissframe$X1,mean)
#' }
#'
#' @keywords survey
#' @export
var.bin <- function(x, bins = 3, iter.max = 100) {
    if (!is(x,"numeric") & !is(x,"integer")) 
        stop("Kmeans clustering not applicable")
    km <- kmeans(x, bins, iter.max = iter.max)
    km$cluster2 <- as.factor(km$cluster)
    levels(km$cluster2) <- levels(km$cluster2)[rank(km$centers)]
    var.binned <- as.numeric(as.character((km$cluster2)))
    return(var.binned)
}
