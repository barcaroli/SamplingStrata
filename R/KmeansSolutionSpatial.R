# ==============================================================================
# File:      KmeansSolutionSpatial.R
# Functions: KmeansSolutionSpatial() [exported]
#
# Purpose:
#   Builds an initial solution for the "spatial" method by clustering the frame
#   units of each domain with k-means on the target variables (on a random
#   subsample of at most maxsamp units), evaluating the strata with
#   aggrStrataSpatial(), i.e. including the spatial autocorrelation component
#   (variogram model of the prediction errors given by vgm_model).
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
#' @param fitting Default is 1, to be used when the frame contains the prediction-error
#'   variances (var columns). A value < 1 (R2 of the prediction model) is meant only for
#'   the rescaling correction, with all variances = 0: with positive variances the function
#'   stops (the error variance would be counted twice). If only R2 is known, use fitting = 1
#'   with a constant variance var(predictions) * (1 - R2) / R2. See buildStrataDFSpatial().
#' @param range Range of the correlation function of the prediction errors
#'   (in the units of the coordinates), as defined by the variogram model given in
#'   vgm_model. 0 means no spatial correlation.
#' @param kappa Factor used in evaluating spatial autocorrelation with the exponential
#'   and gaussian models: 1 if range is the range parameter of the variogram (as
#'   returned by gstat), 3 if it is the practical range. Not used by the spherical model.
#'   Any finite value > 0 is accepted
#'   (a value other than 1 and 3 gives a warning); 0 or negative values give an error.
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
#' @param nugget Relative nugget c0/(c0+c1) of the variogram of the prediction errors,
#'   in [0, 1]: the share of the error variance that is not spatially correlated. Default is 0.
#' @param vgm_model Variogram model fitted on the prediction errors: "exponential"
#'   (default), "spherical" or "gaussian" (or the gstat abbreviations "Exp", "Sph",
#'   "Gau"). A single value, or one value for each target variable Y.
#'   See buildStrataDFSpatial() for the correlation functions.
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
#' # relative nugget c0/(c0+c1)
#' nug <- fit.vgm$var_model$psill[1] / sum(fit.vgm$var_model$psill)
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
#'                             nugget = nug,
#'                             nstrata=NA,
#'                             maxclusters = 5)
#' ############################
#' # Analysis and visualization
#' ############################
#' strataKm <- aggrStrataSpatial(dataset=frame,
#'                         fitting = 1,
#'                         range = fit.vgm$var_model$range[2],
#'                         kappa=1,
#'                         nugget = nug,
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
                                   maxsamp = 5000,
                                   nugget = 0,
                                   vgm_model = "exponential")
{
  cv <- errors
  vgm_model <- vgmModelNames(vgm_model)
  # kappa: any finite positive value is admissible (the correlation is
  # exp(-kappa * d / range), only kappa / range matters); 1 (range parameter of
  # the exponential model) and 3 (practical range) are the usual conventions
  if (length(kappa) != 1 || !is.finite(kappa) || kappa <= 0)
    stop("'kappa' must be a single finite value > 0", call. = FALSE)
  # (kappa is not used by the spherical model)
  if (!kappa %in% c(1, 3) && any(vgm_model != "spherical"))
    warning("kappa = ", kappa, " is unusual: use kappa = 1 when range is the range ",
            "parameter of the variogram (as returned by gstat), kappa = 3 ",
            "when range is the practical range.", call. = FALSE)
  # fitting < 1 together with prediction-error variances counts the error
  # variance twice (see buildStrataDFSpatial)
  vcols <- grep("^var[0-9]+$", colnames(frame), ignore.case = TRUE)
  if (any(fitting < 1) && length(vcols) > 0 &&
      any(unlist(frame[, vcols]) > 0, na.rm = TRUE))
    stop("fitting < 1 cannot be used together with positive prediction-error variances: ",
            "the error variance would be counted twice. Either omit fitting (default 1) ",
            "and keep the prediction-error variances in the frame, or, if only the R2 of ",
            "the prediction model is known, use fitting = 1 with a constant variance ",
            "equal to var(predictions) * (1 - R2) / R2.", call. = FALSE)
  nvariables <- ncol(cv) - 2
  stmt1 <- "solution <- kmeans(stratacorr[,"
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
                              nugget = nugget,
                              vgm_model = vgm_model,
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
    # last point drawn, so that consecutive points are joined by a line
    prev_x <- 1
    prev_y <- sum(v)
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
                                nugget = nugget,
                                vgm_model = vgm_model,
                                vett = solution,
                                dominio = k)
      v <- bethel(aggr, errorscorr, minnumstrat = minnumstrat)
      if (showPlot == TRUE) {
        lines(c(prev_x, i), c(prev_y, sum(v)), type = "b")
        prev_x <- i
        prev_y <- sum(v)
      }
      if (sum(v) <= best[k]) {
        bestsolution <- solution
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
  # build the data frame column by column: cbind() would turn everything into
  # character and as.integer() on the resulting factor would return level codes
  solutionKmean <- data.frame(id = id,
                              suggestions = as.integer(as.character(suggestions)),
                              domainvalue = as.integer(domainvalue),
                              stringsAsFactors = TRUE)
  return(solutionKmean)
}
