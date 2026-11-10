# ==============================================================================
# File:      KmeansSolution2.R
# Functions: KmeansSolution2() [exported]
#
# Purpose:
#   Builds an initial solution for the "continuous" method by clustering the
#   frame units of each domain with k-means on the target variables Y. If
#   nstrata is NA, it also suggests, for each domain, the number of strata (up
#   to maxclusters, by default half the number of units) that minimises the
#   Bethel sample size.
# ==============================================================================

#' @title Initial solution obtained by applying kmeans clustering of frame units
#'
#' @description
#' This function has to be used only in conjunction with "optimizeStrata2",
#' or with "optimStrata" when method = "continuous", i.e. in the
#' case of optimizing with only continuous stratification variables.
#' The function "KmeansSolution2" has a twofold objective:
#' - to give indication about a possible best number of final strata (by fixing a
#' convenient value for "maxclusters", and leaving NA to "nstrata";
#' - to give an initial solution fo the optimization step.
#' If the parameter "nstrata" is not indicated, the optimal number of clusters
#' is determined inside each domain, and the overall solution is obtained by concatenating optimal
#' clusters obtained in domains.
#' The result is a dataframe with two columns: the first indicates the clusters,
#' the second the domains.
#'
#' @param frame The (mandatory) dataframe containing the information related to each unit in the sampling frame.
#' @param model The (optional) dataframe containing the information related to models used to predict
#'   values of the target variables.
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
#' data("swissmunicipalities")
#' swissmunicipalities$id <- c(1:nrow(swissmunicipalities))
#' swissmunicipalities$dom <- 1
#' frame <- buildFrameDF(swissmunicipalities,
#'                       id = "id",
#'                       domainvalue = "REG",
#'                       X = c("Pop020", "Pop2040"),
#'                       Y = c("Pop020", "Pop2040")
#' )
#' cv <- NULL
#' cv$DOM <- "DOM1"
#' cv$CV1 <- 0.1
#' cv$CV2 <- 0.1
#' cv <- as.data.frame(cv)
#' cv <- cv[rep(row.names(cv),7),]
#' cv$domainvalue <- c(1:7)
#' cv
#'
#' # Solution with kmean clustering
#' kmean <- KmeansSolution2(frame,model=NULL,errors=cv,nstrata=NA,maxclusters=4)
#' # number of strata to be obtained in each domain in final solution
#' nstrat <- tapply(kmean$suggestions,
#'                  kmean$domainvalue,
#'                  FUN=function(x) length(unique(x)))
#' # Prepare suggestions for optimization: one chromosome from the k-means
#' # solution, the other from the univariate rule
#' sugg <- prepareSuggestion(frame = frame,
#'                           nStrata = nstrat,
#'                           seeds = c("kmean", "univariate"),
#'                           kmean = kmean)
#' # Optimization
#' # Attention: number of strata must be the same than in kmeans
#' solution <- optimStrata (
#'   method = "continuous",
#'   errors = cv,
#'   framesamp=frame,
#'   iter = 50,
#'   pops = 20,
#'   nStrata = nstrat,
#'   suggestions = sugg,
#'   writeFiles = FALSE,
#'   showPlot = FALSE,
#'   parallel = FALSE
#' )
#' }
#'
#' @keywords survey
#' @export
KmeansSolution2 <- function(frame,
                            model=NULL,
                            errors,
                            nstrata = NA,
                            minnumstrat =2,
                            maxclusters = NA,
                            showPlot = TRUE) {
  colnames(frame) <- toupper(colnames(frame))
  colnames(errors) <- toupper(colnames(errors))
  nvariables <- ncol(errors) - 2
  stmt1 <- "solution <- kmeans(framecorr[,"
  stmt2 <- "c("
  for (i in (1:nvariables)) {
    if (i < nvariables) stmt2 <- paste(stmt2,"'Y",i,"',",sep="")
    if (i == nvariables) stmt2 <- paste(stmt2,"'Y",i,"')",sep="")
  }
  suggestions <- NULL
  domainvalue <- NULL
  id <- NULL
  ndom <- length(unique(frame$DOMAINVALUE))
  solution <- NULL
  best <- NULL
  best_num_strata <- NULL
  for (k in (unique(frame$DOMAINVALUE))) {
    framecorr <- frame[frame$DOMAINVALUE == k,]
    framecorr[,grep("X",colnames(framecorr))] <- NULL
    errorscorr <- errors[errors$DOMAINVALUE == k,]
    stmt3 <- paste("],2)$cluster",sep="")
    stmt <- paste(stmt1,stmt2,stmt3,sep="")
    eval(parse(text=stmt))
    framecorr$X1 <- solution
    aggr <- buildStrataDF(dataset=framecorr,
                           model=NULL,
                           progress=FALSE,
                           verbose=FALSE)
    v <- bethel(aggr, errorscorr, minnumstrat = minnumstrat)

    best[k] <- sum(v)
    if (is.na(maxclusters)) {
      times <- round(nrow(framecorr)*0.5,0)
    }
    if (!is.na(maxclusters)) {
      times <- min(maxclusters,round(nrow(framecorr)*0.5,0))
    }
    if (showPlot == TRUE) {
  	  plot(2,sum(v),xlim=c(1,times),ylim=c(0,1.5*sum(v)),type="p",
         ylab="Sample size",xlab="Number of clusters")
  	  tit <- paste("title('kmeans clustering in domain ",k,"')",sep="")
  	  eval(parse(text=tit))
    }
    bestsolution <- NULL
    if (is.na(nstrata)) {min = 3; max = times }
    if (!is.na(nstrata)) {min = nstrata; max = nstrata }
    for (i in min:max) {
      stmt3 <- paste("],",i,")$cluster",sep="")
      stmt <- paste(stmt1,stmt2,stmt3,sep="")
      eval(parse(text=stmt))
      framecorr$X1 <- solution
      aggr <- buildStrataDF(dataset=framecorr,
                             model=NULL,
                             progress=FALSE,
                             verbose=FALSE)
      v <- bethel(aggr, errorscorr, minnumstrat = minnumstrat)
      # cat("\n",sum(v))
      if (showPlot == TRUE) points(i,sum(v))
      if (sum(v) <= best[k]) {
        bestsolution <- solution
        best_num_strata[k] <- i
        best[k] <- sum(v)
        cat("\nStrata ",i)
        cat("\nSample size ",sum(v))
      }
    }
    bestsolution <- as.factor(bestsolution)
    levels(bestsolution) <- c(1:length(levels(bestsolution)))
    suggestions <- c(suggestions,bestsolution)
    domainvalue <- c(domainvalue,rep(k,nrow(framecorr)))
    id <- c(id,as.character(frame$ID[frame$DOMAINVALUE == k]))
  }
  cat("\n-----------------")
  cat("\n Kmeans solution ")
  cat("\n-----------------")
  for (i in c(1:ndom)) {
    cat("\n *** Domain: ",i," ***")
    cat("\n Number of strata: ",best_num_strata[i])
    cat("\n Sample size     : ",best[i])
  }
  # build the data frame column by column: cbind() would turn everything into
  # character and as.integer() on the resulting factor would return level codes
  solutionKmean <- data.frame(id = id,
                              suggestions = as.integer(as.character(suggestions)),
                              domainvalue = as.integer(domainvalue),
                              stringsAsFactors = TRUE)
  return(solutionKmean)
}
