# ==============================================================================
# File:      KmeansSolutionSpatial.R
# Functions: KmeansSolutionSpatial() [exported]
#
# Purpose:
#   Builds an initial solution for the "spatial" method by clustering the frame
#   units of each domain with k-means on the target variables (on a random
#   subsample of at most maxsamp units), evaluating the strata with
#   aggrStrataSpatial(), i.e. including the spatial autocorrelation component.
#   If nstrata is NA, it also suggests the number of strata that minimises the
#   Bethel sample size.
# ==============================================================================

# Depends on aggrStrata_spatial.r being sourced first: that file is the
# actual SamplingStrata package source for aggrStrataSpatial, fixed in place
# (rejects fitting = 0 with an explanatory error, and fixes the NaN it used
# to produce when range = 0 -- see aggrStrata_spatial.r for the analysis).
#' @title Initial solution obtained by applying kmeans clustering of frame units
#'
#' @description
#' This function has to be used only in conjunction with "optimizeStrataSpatial", i.e. in the
#' case of optimizing with only continuous stratification variables
#' and with a component of spatial autocorrelation.
#' The function "KmeansSolutionSpatial" has a twofold objective:
#' - to give indications about a possible best number of final strata (by fixing a
#' convenient value for "maxclusters", and leaving NA to "nstrata";
#' - to give an initial solution fo the optimization step.
#' If the parameter "nstrata" is not indicated, the optimal number of clusters
#' is determined inside each domain, and the overall solution is obtained by concatenating optimal
#' clusters obtained in domains.
#' The result is a dataframe with two columns: the first indicates the clusters,
#' the second the domains.
#'
#' @param frame The (mandatory) dataframe containing the information related to each unit in the sampling frame.
#' @param fitting Fitting of the model(s).
#'   Default is 1.
#' @param range Maximum range for spatial autocorrelation
#' @param kappa Factor used in evaluating spatial autocorrelation.
#'   Default is 3.
#' @param errors The (mandatory) dataframe containing the precision constraints on target variables.
#' @param nstrata Number of aggregate strata (if NULL, it is optimized by varying the number of cluster
#'   from 2 to half number of atomic strata). Default is NA.
#' @param minnumstrat Minimum number of units to be selected in each stratum. Default is 2.
#' @param maxclusters Maximum number of clusters to be considered in the
#'   execution of kmeans algorithm. If not indicated it will be set equal to the
#'   number of atomic strata divided by 2.
#' @param showPlot Allows to visualise on a plot the different sample sizes for each number of aggregate strata. Default is TRUE.
#' @param maxsamp This is the maximum number of frame units per domain accepted by the function.
#'   If the number of units exceeds this maximum, a subset (size equal to maxsamp) is sampled from the frame.
#'
#' @return
#' A dataframe containing the solution
#'
#' @author
#' Giulio Barcaroli
#'
#' @examples
#' \dontrun{
#' library(sp)
#' library(gstat)
#' library(automap)
#' library(SamplingStrata)
#' #################
#' # data
#' #################
#' # locations (155 observed points)
#' data("meuse")
#' # grid of points (3103)
#' data("meuse.grid")
#' meuse.grid$id <- c(1:nrow(meuse.grid))
#' coordinates(meuse)<-c("x","y")
#' coordinates(meuse.grid)<-c("x","y")
#' #################
#' # kriging
#' #################
#' v <- variogram(lead ~ dist + soil, data=meuse)
#' fit.vgm <- autofitVariogram(lead ~ elev + soil, meuse, model = "Exp")
#' plot(v, fit.vgm$var_model)
#' fit.vgm$var_model
#' # model    psill    range
#' # 1   Nug 1524.895   0.0000
#' # 2   Exp 8275.431 458.3303
#' g <- NULL
#' g <- gstat(g, "Pb", lead ~ dist + soil, meuse)
#' g
#' vm <- variogram(g)
#' vm.fit <- fit.lmc(vm, g, vgm(psill=fit.vgm$var_model$psill[2],
#'                   model="Exp", range=fit.vgm$var_model$range[2],
#'                   nugget=fit.vgm$var_model$psill[1]))
#' # Prediction on the whole grid
#' preds <- predict(vm.fit, meuse.grid)
#' names(preds)
#' # [1] "Pb.pred" "Pb.var"
#' preds$Pb.pred <- ifelse(preds$Pb.pred < 0,0,preds$Pb.pred)
#' df <- NULL
#' df$Pb.pred <- preds@data$Pb.pred
#' df$Pb.var <- preds@data$Pb.var
#' df$dom1 <- 1
#' df <- as.data.frame(df)
#' df$id <- meuse.grid$id
#' #####################################
#' # Optimization with kmeans clustering
#' #####################################
#' frame <- buildFrameDF(df=df,
#'                       id="id",
#'                       X=c("Pb.pred"),
#'                       Y=c("Pb.pred"),
#'                       domainvalue = "dom1")
#' frame$var1 <- df$Pb.var
#' frame$lon <- meuse.grid$x
#' frame$lat <- meuse.grid$y
#' cv <- as.data.frame(list(DOM=rep("DOM1",1),
#'                          CV1=rep(0.05,1),
#'                          domainvalue=c(1:1) ))
#' km <- KmeansSolutionSpatial(frame,
#'                             errors = cv,
#'                             fitting = 1,
#'                             range = fit.vgm$var_model$range[2],
#'                             kappa=1,
#'                             nstrata=NA,
#'                             maxclusters = 5)
#' ############################
#' # Analysis and visualization
#' ############################
#' strataKm <- aggrStrataSpatial(dataset=frame,
#'                         fitting = 1,
#'                         range = fit.vgm$var_model$range[2],
#'                         kappa=1,
#'                         vett=km$suggestions,
#'                         dominio=1)
#' strataKm$SOLUZ <- bethel(strataKm,cv)
#' sum(strataKm$SOLUZ)
#' framenew <- frame
#' framenew$LABEL <- km$suggestions
#' strataKm$STRATO <- strataKm$stratum
#' ssKm <- summaryStrata(framenew,strataKm)
#' ssKm
#' frameres <- SpatialPointsDataFrame(data=framenew,
#'                                    coords=cbind(framenew$lon,framenew$lat) )
#' frameres2 <- SpatialPixelsDataFrame(points=frameres[c("lon","lat")],
#'                                     data=framenew)
#' frameres2$LABEL <- as.factor(frameres2$LABEL)
#' spplot(frameres2,c("LABEL"), col.regions=bpy.colors(5))
#' }
#'
#' @keywords survey
#' @export
KmeansSolutionSpatial <- function (frame,
                                   fitting = 1, 
                                   range = c(0), 
                                   kappa = 3,
                                   errors, 
                                   nstrata = NA, 
                                   minnumstrat = 2,
                                   maxclusters = NA,
                                   showPlot = TRUE,
                                   maxsamp = 5000)
{
  cv <- errors
  nvariables <- ncol(cv) - 2
  stmt1 <- "solution <- kmeans(frame[,"
  stmt2 <- "c("
  for (i in (1:nvariables)) {
    if (i < nvariables) 
      stmt2 <- paste(stmt2, "'Y", i, "',", sep = "")
    if (i == nvariables) 
      stmt2 <- paste(stmt2, "'Y", i, "')", sep = "")
  }
  suggestions <- NULL
  domainvalue <- NULL
  id <- NULL
  ndom <- nrow(cv)
  solution <- NULL
  best <- rep(0, ndom)
  best_num_strata <- rep(0, ndom)
  for (k in (1:nrow(cv))) {
    # Subsample large domains before the spatial computation: aggrStrataSpatial
    # builds full n x n pairwise distance/covariance matrices (O(n^2) memory),
    # so a domain with tens of thousands of points exhausts memory ("cannot
    # allocate vector"). The K-means step only provides an initial number of
    # strata, so a random subsample of up to 'maxsamp' points is sufficient.
    # Set maxsamp = Inf to disable. Below the cap, behaviour is unchanged.
    idx_k <- which(frame$domainvalue == k)
    if (length(idx_k) > maxsamp) idx_k <- sample(idx_k, maxsamp)
    stratacorr <- frame[idx_k, ]
    errorscorr <- errors[errors$domainvalue == k, ]
    aggr <- aggrStrataSpatial(dataset = stratacorr,
                              fitting = fitting,
                              range = range,
                              kappa = kappa,
                              vett = rep(1, nrow(stratacorr)),
                              dominio = k)
    v <- bethel(aggr, errorscorr, minnumstrat = minnumstrat)
    sum(v)
    best[k] <- sum(v)
    if (is.na(maxclusters)) {
      times <- round(nrow(stratacorr) * 0.5, 0)
    }
    if (!is.na(maxclusters)) {
      times <- min(maxclusters, round(nrow(stratacorr) * 
                                        0.5, 0))
    }
    if (showPlot == TRUE) {
      plot(1, sum(v), xlim = c(1, times), ylim = c(0, 
                                                   1.5 * sum(v)), type = "p", ylab = "Sample size", 
           xlab = "Number of clusters")
      tit <- paste("title('kmeans clustering in domain ", 
                   k, "')", sep = "")
      eval(parse(text = tit))
    }
    bestsolution <- NULL
    if (is.na(nstrata)) {
      min = 2
      max = times
    }
    if (!is.na(nstrata)) {
      min = nstrata
      max = nstrata
    }
    for (i in min:max) {
      stmt3 <- paste("],", i, ")$cluster", sep = "")
      stmt <- paste(stmt1, stmt2, stmt3, sep = "")
      eval(parse(text = stmt))
      aggr <- aggrStrataSpatial(dataset = stratacorr, 
                                fitting = fitting, 
                                range = range, 
                                kappa = kappa,
                                vett = solution[idx_k],
                                dominio = k)
      v <- bethel(aggr, errorscorr, minnumstrat = minnumstrat)
      if (showPlot == TRUE) 
        points(i, sum(v))
      if (sum(v) <= best[k]) {
        bestsolution <- solution[idx_k]
        best_num_strata[k] <- i
        best[k] <- sum(v)
      }
    }
    bestsolution <- as.factor(bestsolution)
    levels(bestsolution) <- c(1:length(levels(bestsolution)))
    suggestions <- c(suggestions, bestsolution)
    domainvalue <- c(domainvalue, rep(k, nrow(stratacorr)))
    id <- c(id,as.character(frame$id[idx_k]))
  }
  cat("\n-----------------")
  cat("\n Kmeans solution ")
  cat("\n-----------------")
  for (i in c(1:ndom)) {
    cat("\n *** Domain: ", i, " ***")
    cat("\n Number of strata: ", best_num_strata[i])
    cat("\n Sample size     : ", best[i])
  }
  solutionKmean <- as.data.frame(cbind(id,suggestions,domainvalue),stringsAsFactors = TRUE)
  solutionKmean$domainvalue <- as.integer(solutionKmean$domainvalue)
  return(solutionKmean)
}
