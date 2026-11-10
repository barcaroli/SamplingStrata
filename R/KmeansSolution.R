# ==============================================================================
# File:      KmeansSolution.R
# Functions: KmeansSolution() [exported]
#
# Purpose:
#   Builds an initial solution ("suggestion") for optimizeStrata() (method
#   "atomic") by clustering the atomic strata with k-means on the means of the
#   target variables. If the number of strata is not given, it is chosen, for
#   each domain, as the value (up to maxclusters) that minimises the Bethel
#   sample size.
# ==============================================================================

#' @title Initial solution obtained by applying kmeans clustering of atomic strata
#'
#' @description
#' In order to speed the convergence towards the optimal solution, an initial one can be
#' given as "suggestion" to "optimizeStrata" function.
#' The function "KmeansSolution" produces this initial solution using the k-means
#' algorithm by clustering atomic
#' strata on the basis of the values of the means of target variables in them.
#' Also, if the parameter "nstrata" is not indicated, the optimal number of clusters
#' is determined inside each domain, and the overall solution is obtained by concatenating optimal
#' clusters obtained in domains.
#' The result is a dataframe with two columns: the first indicates the clusters,
#' the second the domains.
#'
#' @param strata The (mandatory) dataframe containing the information related to atomic strata.
#' @param errors The (mandatory) dataframe containing the precision constraints on target variables.
#' @param nstrata Number of aggregate strata (if NULL, it is optimized by varying the number of cluster
#'   from 2 to half number of atomic strata). Default is NA.
#' @param minnumstrat Minimum number of units to be selected in each stratum. Default is 2.
#' @param maxclusters Maximum number of clusters to be considered in the
#'   execution of kmeans algorithm. If not indicated it will be set equal to the
#'   number of atomic strata divided by 2.
#' @param showPlot Allows to visualise on a plot the different sample sizes for each number of aggregate strata. Default is TRUE.
#'
#' @return
#' A dataframe containing the solution
#'
#' @author
#' Giulio Barcaroli
#'
#' @examples
#' \dontrun{
#' library(SamplingStrata)
#' data(swisserrors)
#' data(swissframe)
#' # atomic strata of the frame (the same that optimStrata() builds internally)
#' strata <- buildStrataDF(swissframe, progress = FALSE)
#' # suggestion
#' solutionKmean <- KmeansSolution(strata = strata,
#'                 errors = swisserrors,
#'                 nstrata = NA,
#'                 maxclusters = 10,
#'                 showPlot = TRUE)
#'
#' # number of strata to be obtained in each domain
#' nstrat <- tapply(solutionKmean$suggestions,
#'                  solutionKmean$domainvalue,
#'                  FUN=function(x) length(unique(x)))
#'
#' # optimisation of sampling strata
#' solution <- optimStrata (
#'   method = "atomic",
#'   errors = swisserrors,
#'   framesamp = swissframe,
#'   nStrata = nstrat,
#'   suggestions = solutionKmean,
#'   iter = 50,
#'   pops = 10,
#'   parallel = FALSE
#' )
#' }
#'
#' @keywords survey
#' @export
KmeansSolution <- function(strata,
                             errors,
                             nstrata = NA,
                             minnumstrat =2,
                             maxclusters = NA,
                             showPlot = TRUE) {
  nvariables <- ncol(errors) - 2
  stmt1 <- "solution <- kmeans(stratacorr[,"
  stmt2 <- "c("
  for (i in (1:nvariables)) {
    if (i < nvariables) stmt2 <- paste(stmt2,"'M",i,"',",sep="")
    if (i == nvariables) stmt2 <- paste(stmt2,"'M",i,"')",sep="")
  }
  suggestions <- NULL
  domainvalue <- NULL
  ndom <- length(unique(strata$DOM1))
  solution <- NULL
  best <- rep(0,ndom)
  best_num_strata <- rep(0,ndom)
    for (k in (unique(strata$DOM1))) {
      stratacorr <- strata[strata$DOM1 == k,]
      errorscorr <- errors[errors$domainvalue == k,]
      aggr <- aggrStrata(strata=stratacorr,
                         nvar=nvariables,
                         censiti=0,
                         vett=rep(1,nrow(stratacorr)),
                         dominio=k)
      v <- bethel(aggr, errorscorr, minnumstrat = minnumstrat)
      sum(v)
      best[k] <- sum(v)
      if (is.na(maxclusters)) {
        times <- round(nrow(stratacorr)*0.5,0)
      }
      if (!is.na(maxclusters)) {
        times <- min(maxclusters,round(nrow(stratacorr)*0.5,0))
      }
      if (showPlot == TRUE) {
    	  plot(1,sum(v),xlim=c(1,times),ylim=c(0,1.5*sum(v)),type="p",
           ylab="Sample size",xlab="Number of clusters")
    	  tit <- paste("title('kmeans clustering in domain ",k,"')",sep="")
    	  eval(parse(text=tit))
      }
      bestsolution <- NULL
      if (is.na(nstrata)) {min = 2; max = times }
      if (!is.na(nstrata)) {min = nstrata; max = nstrata }
      for (i in min:max) {
        stmt3 <- paste("],",i,")$cluster",sep="")
        stmt <- paste(stmt1,stmt2,stmt3,sep="")
        solution <- tryCatch(eval(parse(text = stmt)))
        aggr <- aggrStrata(strata=stratacorr,
                           nvar=nvariables,
                           censiti=0,
                           vett=solution,
                           dominio=k)
        v <- bethel(aggr, errorscorr, minnumstrat = minnumstrat)
        # cat("\n",sum(v))
        if (showPlot == TRUE) points(i,sum(v))
        if (sum(v) <= best[k]) {
          bestsolution <- solution
          best_num_strata[k] <- i
          best[k] <- sum(v)
        }
      }
      suggestions <- c(suggestions,bestsolution)
      domainvalue <- c(domainvalue,rep(k,nrow(stratacorr)))
  }
  solutionKmean <- as.data.frame(cbind(suggestions,domainvalue),stringsAsFactors = TRUE)
  solutionKmean$domainvalue <- as.integer(solutionKmean$domainvalue)
  totsize <- 0
  cat("\n-------------------")
  cat("\n  Kmeans solution ")
  cat("\n-------------------")
  for (i in (unique(strata$DOM1))) {
    cat("\n *** Domain: ", i, " ***")
    cat("\n Number of strata: ", best_num_strata[i])
    cat("\n Sample size     : ", best[i])
    totsize <- totsize + best[i]
  }
  cat("\n-------------------")
  cat("\n Total size: ",totsize)
  cat("\n-------------------")
  solutionKmean
}