# ==============================================================================
# File:      assignStrataLabel.R
# Functions: assignStrataLabel() [exported]
#
# Purpose:
#   Assigns the optimized stratum labels to new units of a frame (not used in
#   the optimisation), using the strata structure (ranges of the X variables)
#   returned by summaryStrata() after an optimisation with optimizeStrata2() /
#   optimStrata(method = "continuous").
# ==============================================================================

# --------------------------------------------------
# assignStrataLabel
# Function to assign the optimized strata labels
# to new sampling units in the frame
# 18 April 2019
# Giulio Barcaroli
# --------------------------------------------------

#' @title Function to assign the optimized strata labels
#'
#' @description
#' Function to assign the optimized strata labels to new sampling units in the frame
#' on the basis of the strata structure obtained by executing the function
#' 'summaryStrata' after optimizing with 'optimizeStrata2'
#'
#' @param dataset dataset with new sampling units in the frame
#' @param s structure of the strata
#'
#' @return
#' The same dataset in input with the label of the optimized stratum
#'
#' @examples
#' \dontrun{
#' library(SamplingStrata)
#' data("swissmunicipalities")
#' data("errors")
#' errors$CV1 <- 0.1
#' errors$CV2 <- 0.1
#' errors <- errors[rep(row.names(errors),7),]
#' errors$domainvalue <- c(1:7)
#' errors
#' swissmunicipalities$id <- c(1:nrow(swissmunicipalities))
#' swissmunicipalities$domain = 1
#' frame <- buildFrameDF(swissmunicipalities,
#'                       id = "id",
#'                       domainvalue = "REG",
#'                       X = c("Surfacesbois","Surfacescult"),
#'                       Y = c("Pop020", "Pop2040")
#' )
#' solution <- optimizeStrata2 (
#'   errors,
#'   frame,
#'   nStrata = 5,
#'   iter = 10,
#'   pops = 10,
#'   writeFiles = FALSE,
#'   showPlot = TRUE,
#'   parallel = FALSE)
#' strataStructure <- summaryStrata(solution$framenew, solution$aggr_strata)
#' strataStructure
#'
#' newset <- assignStrataLabel(solution$framenew,strataStructure)
#' }
#'
#' @keywords survey
#' @export
assignStrataLabel <- function(dataset, s) 
  {
  colnames(dataset) <- toupper(colnames(dataset))
  colnames(s) <- toupper(colnames(s))
  # check on dataset content
  if (!"DOMAINVALUE" %in% colnames(dataset)) stop("No DOMAINVALUE variable in dataset")
  nvarX_1 <- length(grep("LOWER",colnames(s)))
  nvarX_2 <- length(grep("X",colnames(dataset))) 
  if (nvarX_1 != nvarX_2)  stop("Number of X variables in dataset not compatible with strata structure")
  df <- NULL  
  for (dom in unique(dataset$DOMAINVALUE)) {
    for (strat in unique(s$STRATUM[s$DOMAIN==dom])) {
      cat("\nDomain ",dom,"  Stratum ",strat)
      stmt <- "strato <- dataset[dataset$DOMAINVALUE == dom & ("
      # conditions for inclusion
      for (i in (1:nvarX_1)) {
        if (i != nvarX_1) {
          stmt <- paste(stmt,"(dataset$X",i," >= s[s$DOMAIN==dom & s$STRATUM==strat,'LOWER_X",i,"'] & dataset$X",i," <= s[s$DOMAIN==dom & s$STRATUM==strat,'UPPER_X",i,"']) | ",sep="")
      }                    
        if (i == nvarX_1) {
          stmt <- paste(stmt,"(dataset$X",i," >= s[s$DOMAIN==dom & s$STRATUM==strat,'LOWER_X",i,"'] & dataset$X",i," <= s[s$DOMAIN==dom & s$STRATUM==strat,'UPPER_X",i,"']) )",sep="")
        }
      }
      # conditions for exclusion
      stmt <- paste(stmt," & (")
      for (i in (1:nvarX_1)) {
        if (i != nvarX_1) {
          stmt <- paste(stmt,"!(dataset$X",i," < s[s$DOMAIN==dom & s$STRATUM==strat,'LOWER_X",i,"']) & !(dataset$X",i," > s[s$DOMAIN==dom & s$STRATUM==strat,'UPPER_X1']) & ",sep="")
        }                    
        if (i == nvarX_1) {
          stmt <- paste(stmt,"!(dataset$X",i," < s[s$DOMAIN==dom & s$STRATUM==strat,'LOWER_X",i,"']) & !(dataset$X",i," > s[s$DOMAIN==dom & s$STRATUM==strat,'UPPER_X",i,"']) ), ]",sep="")
        }
      }
      eval(parse(text=stmt))
      strato$STRATO <- paste(dom,strat,sep="***")
      stratomin <- strato[!(strato$ID %in% df$ID),]
      df <- rbind(df,stratomin)
    }
  }
  # t <- as.data.frame(table(df$STRATO))
  # write.table(t,"strati_out.csv",sep=";",row.names = F)

  return(df)
}

