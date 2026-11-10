# ==============================================================================
# File:      buildFrameDF.R
# Functions: buildFrameDF() [exported]
#
# Purpose:
#   Builds the "frame" dataframe in the standard format required by the
#   package from a generic dataset: renames the identifier to id, the
#   stratification variables to X1, X2, ..., the target variables to Y1, Y2, ...
#   and the domain variable to domainvalue (converted to integer). If a weight
#   variable is given (frame derived from a sample survey), each unit is
#   replicated a number of times equal to its rounded weight.
# ==============================================================================

#' @title Builds the "sampling frame" dataframe from a dataset containing information on all the units in the population of reference
#'
#' @description
#' This function allows to build the information regarding the sampling frame of the population of reference.
#' Mandatory variables are:
#' (i) the name of the dataset containing the sampling frame of the population of reference
#' (ii) an identifier (Id)
#' (iii) a set of auxiliary variables X
#' (iv) a set of target variables Y
#' (v) the indicator of the domain to which the unit belongs
#' A further optional parameter is the the name of the variable containing weights (in case instead of a frame we use a sampling survey)
#'
#' @param df This is the name of the dataframe containing the information
#'   on all the units in population of reference.
#' @param id This is the name of the identifier in the sampling frame.
#' @param X A character vector containing the names of the auxiliary variables
#'   in the frame dataset
#' @param Y A character vector containing the names of the target variables
#'   in the frame dataset
#' @param domainvalue The name of the variable in the frame dataset that contains the
#'   indication of the domains to which the units belong.
#' @param weight The name of the variable in the frame dataset that contains the
#'   weights (in case we do not use a proper frame, but a sampling survey.
#'
#' @return
#' A dataframe
#'
#' @author
#' Giulio Barcaroli
#'
#' @examples
#' \dontrun{
#' data(swissmunicipalities)
#' id = "Nom"
#' X = c("Surfacesbois","Surfacescult")
#' Y = c("Pop020","Pop2040")
#' domainvalue = "REG"
#' frame <- buildFrameDF(swissmunicipalities,id,X,Y,domainvalue)
#' head(frame)
#' }
#'
#' @keywords survey
#' @export
buildFrameDF <- function (df, id, X, Y, domainvalue, weight = NULL) 
{
  vars <- colnames(df)
  if (!id %in% vars) 
    stop("Id name not in the frame variables")
  for (i in (1:length(X))) {
    var = X[i]
    if (!var %in% vars) {
      msg = paste("X variable ", var, " not in the frame variables", 
                  sep = "")
      stop(msg)
    }
  }
  for (i in (1:length(Y))) {
    var = Y[i]
    if (!var %in% vars) {
      msg = paste("Y variable ", var, " not in the frame variables", 
                  sep = "")
      stop(msg)
    }
  }
  if (!domainvalue %in% vars) 
    stop("Domain name not in the frame variables")
  dframe <- NULL
  stmt <- paste("dframe$id <- df$", as.character(id), sep = "")
  eval(parse(text = stmt))
  for (i in (1:length(X))) {
    aux <- paste("dframe$X", i, " <- df$", as.character(X[i]), 
                 sep = "")
    eval(parse(text = aux))
  }
  for (i in (1:length(Y))) {
    target <- paste("dframe$Y", i, " <- df$", as.character(Y[i]), 
                    sep = "")
    eval(parse(text = target))
  }
  stmt <- paste("dframe$domainvalue <- df$", as.character(domainvalue), 
                sep = "")
  eval(parse(text = stmt))
  dframe <- as.data.frame(dframe, stringsAsFactors = TRUE)
  dframe$domainvalue <- as.integer(dframe$domainvalue)
  if (!is.null(weight)) {
    stmt <- paste("dframe$WEIGHT <- df$", as.character(weight), 
                  sep = "")
    eval(parse(text = stmt))
    dframe <- dframe[rep(1:nrow(dframe), round(dframe$WEIGHT)),]
    dframe$WEIGHT <- NULL
  }
  return(dframe)
}