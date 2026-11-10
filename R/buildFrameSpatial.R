# ==============================================================================
# File:      buildFrameSpatial.R
# Functions: buildFrameSpatial() [exported]
#
# Purpose:
#   Builds the "frame" dataframe for the "spatial" method: as buildFrameDF(),
#   plus the prediction variances of the target variables (columns var1, var2,
#   ...) and the geographical coordinates (lon, lat) of each unit.
# ==============================================================================

#' @title Builds the "sampling frame" dataframe from a dataset containing information all the units in the population of reference including spatial
#'
#' @description
#' This function allows to build the information regarding the sampling frame of the population of reference.
#' Mandatory variables are:
#' (i) the name of the dataset containing the sampling frame of the population of reference
#' (ii) an identifier (Id)
#' (iii) a set of auxiliary variables 'X'
#' (iv) a set of target variables 'Y'
#' (v) a set of prediction errors variables 'variance'
#' (vi) longitude
#' (vii) latitude
#' (viii) the indicator of the domain to which the unit belongs
#'
#' @param df This is the name of the dataframe containing the information
#'   on all the units in population of reference.
#' @param id This is the name of the identifier in the sampling frame.
#' @param X A character vector containing the names of the auxiliary variables
#'   in the frame dataset
#' @param Y A character vector containing the names of the target variables
#'   in the frame dataset
#' @param variance A character vector containing the names of the prediction error variables
#'   in the frame dataset
#' @param lon Longitude of the unit
#' @param lat Latitude of the unit
#' @param domainvalue The name of the variable in the frame dataset that contains the
#'   indication of the domains to which the units belong.
#'
#' @return
#' A dataframe
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
#' data("meuse")
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
#' g <- NULL
#' g <- gstat(g, "Pb", lead ~ dist + soil, meuse)
#' g
#' vm <- variogram(g)
#' vm.fit <- fit.lmc(vm, g, vgm(psill=fit.vgm$var_model$psill[2],
#'                              model="Exp", range=fit.vgm$var_model$range[2],
#'                              nugget=fit.vgm$var_model$psill[1]))
#' # Prediction on the whole grid
#' preds <- predict(vm.fit, meuse.grid)
#' names(preds)
#' # [1] "Pb.pred" "Pb.var"
#' preds$Pb.pred <- ifelse(preds$Pb.pred < 0,0,preds$Pb.pred)
#' df <- NULL
#' df$id <- meuse.grid$id
#' df$Pb.pred <- preds@data$Pb.pred
#' df$Pb.var <- preds@data$Pb.var
#' df$lon <- meuse.grid$x
#' df$lat <- meuse.grid$y
#' df$dom1 <- 1
#' df <- as.data.frame(df)
#' frame <- buildFrameSpatial(df=df,
#'                       id="id",
#'                       X=c("Pb.pred"),
#'                       Y=c("Pb.pred"),
#'                       variance=c("Pb.var"),
#'                       lon="lon",
#'                       lat="lat",
#'                       domainvalue = "dom1")
#' head(frame)
#' }
#'
#' @keywords survey
#' @export
buildFrameSpatial <- function(df,id,X,Y,variance,lon,lat,domainvalue) {
  vars <- colnames(df)
  if (!id %in% vars) stop("Id name not in the frame variables")
  if (!lon %in% vars) stop("lon name not in the frame variables")
  if (!lat %in% vars) stop("lat name not in the frame variables")
  
  for (i in (1:length(X))) {
    var = X[i]
    if (!var %in% vars) {
      msg = paste("X variable ",var," not in the frame variables",sep="")
      stop(msg)
    }
  }
  for (i in (1:length(Y))) {
    var = Y[i]
    if (!var %in% vars) {
       msg = paste("Y variable ",var," not in the frame variables",sep="")
       stop(msg)
    }
  }
  if (length(variance) != length(Y)) stop("Number of variances not equal to number of Y's")
  for (i in (1:length(variance))) {
    var = variance[i]
    if (!var %in% vars) {
      msg = paste("Variance variable ",var," not in the frame variables",sep="")
      stop(msg)
    }
  }
  if (!domainvalue %in% vars) stop("Domain name not in the frame variables")  
  dframe <- NULL
  stmt <- paste("dframe$id <- df$",as.character(id),sep="")
  eval(parse(text=stmt))
  for (i in (1:length(X))) {
    aux <- paste("dframe$X",i," <- df$",as.character(X[i]),sep="")
    eval(parse(text=aux))
  }
  for (i in (1:length(Y))) {
    target <- paste("dframe$Y",i," <- df$",as.character(Y[i]),sep="")
    eval(parse(text=target))
  }
  for (i in (1:length(variance))) {
    target <- paste("dframe$var",i," <- df$",as.character(variance[i]),sep="")
    eval(parse(text=target))
  }
  target <- paste("dframe$lon <- df$",as.character(lon),sep="")
  eval(parse(text=target))
  target <- paste("dframe$lat <- df$",as.character(lat),sep="")
  eval(parse(text=target))
  stmt <- paste("dframe$domainvalue <- df$",as.character(domainvalue),sep="")
  eval(parse(text=stmt))
  dframe <- as.data.frame(dframe,stringsAsFactors = TRUE)
  dframe$domainvalue <- as.integer(dframe$domainvalue)
  return(dframe)
}