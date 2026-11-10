# ==============================================================================
# File:      expected_CV.R
# Functions: expected_CV() [exported]
#
# Purpose:
#   Computes the expected coefficients of variation of the target variables in
#   each domain for the aggregated strata of an optimized solution (using the
#   allocation in column SOLUZ), to be compared with the precision
#   constraints.
# ==============================================================================

#' @title Expected coefficients of variation of target variables Y
#'
#' @description
#' Once optimized the sampling frame, this function allows to calculate the expected coefficients of variation on the aggregated strata of the current optimized solution, and to compare them to the precision constraints.
#'
#' @param strata Aggregated strata in the solution obtained by the execution of the 'optimizedStrata' function
#'
#' @return
#' Matrix containing values of the CV's in the different domains
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
#' # calculate CV's on Y's
#' expected_CV(solution$aggr_strata)
#' # compare to precision constraints
#' swisserrors
#' }
#'
#' @keywords survey
#' @export
expected_CV <- function (strata) {
  M_h <- S_h <- NULL
  if (is.null(strata$SOLUZ)) stop("There is no allocation of units in strata")
  ndom <- length(unique(strata$DOM1))
  nvarX <- length(grep("X",colnames(strata)))
  nstr <- length(grep("STR",colnames(strata)))
  nvars <- (ncol(strata) - nstr - 5 - nvarX) / 2
  cv <- matrix(NA,nrow=ndom,ncol=nvars)
  colnames(cv) <- paste("cv(Y",c(1:nvars),")",sep="")
  rownames(cv) <- paste("DOM",c(1:ndom),sep="")
  k<-0
  for (i in (as.numeric(levels(as.factor(strata$DOM1))))) {
    k<-k+1
    stratadom <- strata[strata$DOM1 == i,]
    for (j in 1:nvars) {
      n_h <- stratadom$SOLUZ
      N_h <- stratadom$N
      stmt <- paste("S_h <- stratadom$S",j,sep="")
      eval(parse(text=stmt))
      stmt <- paste("M_h <- stratadom$M",j,sep="")
      eval(parse(text=stmt))
      Y_h <- N_h * M_h
      Var_h <- (N_h^2) * (1 - n_h/N_h) * ((S_h^2)/n_h)
      CV <- sqrt(sum(Var_h)) / sum(Y_h)
      cv[k,j] <- CV 
    }
  }
  cv <- round(cv,7)
  cv
}