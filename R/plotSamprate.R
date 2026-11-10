# ==============================================================================
# File:      plotSamprate.R
# Functions: plotSamprate() [exported]
#
# Purpose:
#   Draws, for each domain, a bar plot of the sampling rates (allocated units /
#   population) in the optimized strata.
# ==============================================================================

#' @title Plotting sampling rates in the different strata for each domain in the solution.
#'
#' @description
#' Once the optimization step has been carried out, by applying this function it is
#' possible to obtain the visualization of the proportion of sampling units in the
#' different strata for each domain in the obtained solution.
#'
#' @param solution Solution obtained by executing optimizeStrata
#' @param dom Identification of the domain
#'
#' @return
#' Plot
#'
#' @examples
#' \dontrun{
#' library(SamplingStrata)
#' data(swisserrors)
#' data(swissstrata)
#' # optimisation of sampling strata
#' solution <- optimizeStrata (
#'     errors = swisserrors,
#'     strata = swissstrata,
#' )
#' # plot of the sampling rates in strata
#' for (i in (1:length(unique(swissstrata$DOM1)))) plotSamprate(solution, i)
#' }
#'
#' @keywords survey
#' @export
plotSamprate <- function (solution, dom) {
  df <- as.data.frame(solution$aggr_strata[solution$aggr_strata$DOM1 == dom,])
#   plot(df$SOLUZ/df$N,
#        ylim = c(0,1),
#        ylab="sampling units / total units",
#        xlab="strata",
# 	   type="h")
  barplot(SOLUZ/N ~ STRATO,data=df, ylab="sampling units / total units",xlab="strata",)
  title(paste("Sampling rate per stratum in domain ",dom,sep=""))
}    