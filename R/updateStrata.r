# ==============================================================================
# File:      updateStrata.r
# Functions: updateStrata() [exported]
#
# Purpose:
#   Assigns to each atomic stratum the label of the optimized aggregated stratum
#   it belongs to ("newstrata"), and optionally writes the result and the
#   structure of the aggregated strata to files.
# ==============================================================================

# ----------------------------------------------------
# Function for assigning new labels to initial strata
# and to report the structure of resulting
# aggregated strata
# Author: Giulio Barcaroli
# Date: 4 January 2012
# Last update: 12 January 2023
# ----------------------------------------------------
#' @title Assigns new labels to atomic strata on the basis of the optimized aggregated strata
#'
#' @description
#' Once optimal stratification has been obtained ('outstrata'), then we need to attribute
#' new strata labels to each atomic stratum.
#' By executing this function, a new dataframe "newstrata" will be obtained with the same structure of the old,
#' ("strata") but with the addition of a new stratum label.
#' By indicating "YES" to "writeFile" parameter, the dataframe "newstrata" will be written to a delimited file
#' ("newstrata.txt").
#' Also a second delimited file ("strata_aggregation.txt") will be outputted,
#' containing the indication of the relations bewteen atomic and aggregated strata.
#'
#' @param strata This is the (mandatory) dataframe containing the information related to the atomic strata to which the
#'   optimisation has been applied to.
#' @param solution List obtained by the execution of the "optimizeStrata" function. The first element of the list is the vector
#'   of the indices corresponding to the optimal solution.
#' @param writeFiles Indicates if at the end of the processing the resulting strata will be outputted in a delimited file.
#'   Default is "FALSE".
#' @param outputFolder A path to a directory ot be used for the outputs. Default is "output"
#'
#' @return
#' A dataframe containing the strata
#'
#' @author
#' Giulio Barcaroli
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
#' # updating sampling strata with new strata labels
#' newstrata <- updateStrata(swissstrata, solution, writeFiles = TRUE)
#' }
#'
#' @keywords survey
#' @export
updateStrata <- function (strata, 
                          solution, 
                          writeFiles = FALSE,
                          outputFolder = file.path(getwd(),"output"))
{
  if (writeFiles == TRUE) {
    # if(dir.exists(outputFolder)){
    #   warning("Folder ", outputFolder," already existed and has been deleted.")
    #   unlink(outputFolder)
    # } 
    if(!dir.exists(outputFolder)) dir.create(outputFolder)
  }
    colnames(strata) <- toupper(colnames(strata))
    newstrata <- strata
    newstrata$AGGR_STRATUM <- solution[[1]]
    ndom <- length(levels(as.factor(strata$DOM1)))
    nvarX <- length(grep("X", names(strata)))
    matstrata <- NULL
    stmt <- "matstrata <- as.data.frame(cbind(newstrata$DOM1,newstrata$AGGR_STRATUM,"
    stmt2 <- "colnames(matstrata) <- c('DOM1','AGGR_STRATUM',"
    stmt3 <- NULL
    if (nvarX > 1) {
        for (i in 1:(nvarX - 1)) {
            stmt <- paste(stmt, "newstrata$X", i, ",", sep = "")
            stmt2 <- paste(stmt2, "'X", i, "',", sep = "")
            stmt3 <- paste(stmt3, "matstrata$X", i, ",", sep = "")
        }
        stmt <- paste(stmt, "newstrata$X", nvarX, "), stringsAsFactors = TRUE)", sep = "")
        eval(parse(text = stmt))
        stmt2 <- paste(stmt2, "'X", nvarX, "')", sep = "")
        eval(parse(text = stmt2))
        stmt3 <- paste(stmt3, "matstrata$X", nvarX, sep = "")
        statement <- paste("matstrord <- matstrata[order(matstrata$DOM1,matstrata$AGGR_STRATUM,", 
            stmt3, "),]", sep = "")
        eval(parse(text = statement))
    }
    if (nvarX == 1) {
		matstrata <- as.data.frame(cbind(newstrata$DOM1,newstrata$AGGR_STRATUM,newstrata$X1),stringsAsFactors = TRUE)
		colnames(matstrata) <- c('DOM1','AGGR_STRATUM','X1')
        matstrord <- matstrata[order(matstrata$DOM1, matstrata$AGGR_STRATUM, 
            matstrata$X1), ]
    }
    if (nvarX == 1) 
        newstrata$STRATUM <- newstrata$X1
    if (nvarX > 1) {
        stmt <- NULL
        stmt <- "newstrata$STRATUM <- paste("
        for (i in 1:(nvarX - 1)) {
            if (i > 0) 
                stmt <- paste(stmt, "newstrata$X", i, ",", sep = "")
        }
        stmt <- paste(stmt, "newstrata$X", nvarX, ",sep='*')", 
            sep = "")
        eval(parse(text = stmt))
    }
    colnames(newstrata)[ncol(newstrata) - 1] <- c("LABEL")
    colnames(newstrata) <- toupper(colnames(newstrata))
    if (writeFiles == TRUE) 
        write.table(newstrata, file = file.path(outputFolder,"newstrata.txt"), sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote = FALSE)
    if (writeFiles == TRUE) 
        write.table(matstrord, file = file.path(outputFolder,"strata_aggregation.txt"), 
            sep = "\t", row.names = FALSE, col.names = TRUE, 
            quote = FALSE)
    # if (writeFiles == TRUE) {
    #   setwd(dire)
    # }
    return(newstrata)
}
