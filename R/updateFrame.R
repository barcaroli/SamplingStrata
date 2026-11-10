# ==============================================================================
# File:      updateFrame.R
# Functions: updateFrame() [exported]
#
# Purpose:
#   Adds to the sampling frame the label of the optimized stratum of each unit,
#   taken from the atomic strata relabelled by updateStrata().
# ==============================================================================

# ----------------------------------------------------
# Function to assign labels of optimal aggregated strata
# to sampling frame
# Author: Giulio Barcaroli
# Date: 4 January 2012
# Last update: 12 January 2023
# ----------------------------------------------------
#' @title Updates the initial frame on the basis of the optimized stratification
#'
#' @description
#' Once optimal stratification has been obtained, and new labels have been attributed to initial
#' atomic strata ("newstrata"), it is important to report the new classification of units in the
#' sampling frame by attributing new strata labels to each unit.
#' By executing this function, a new frame will be obtained with the same structure of the old,
#' but with the addition of a new stratum label.
#' The initial frame must contain a variable named 'domainvalue' that indicates the
#' same values of the domain that has been used with the 'optimizeStrata' function.
#' If no domains have been defined, this variable will contains all 1's, but it must exist
#'
#' @param frame This is the (mandatory) dataframe containing the sampling frame.
#' @param newstrata This is the (mandatory) dataframe containing the information related to the optimisation applied to
#'   initial stratification (new labels applied to atomic strata). It is produced by executing
#'   the "updateStrata" function.
#' @param writeFiles Flag to write or not the new sampling frame into the working directory.
#'   Default is "FALSE"
#' @param outputFolder A path to a directory ot be used for the outputs. Default is "output"
#'
#' @return
#' A dataframe containing the frame
#'
#' @author
#' Giulio Barcaroli
#'
#' @examples
#' #
#' # The following example is realistic, but is time consuming
#' #
#' \dontrun{
#' library(SamplingStrata)
#' data(swisserrors)
#' data(swissstrata)
#' # optimisation of sampling strata
#' solution <- optimizeStrata (
#'     errors = swisserrors,
#'     strata = swissstrata)
#' # updating sampling strata with new strata labels
#' newstrata <- updateStrata(swissstrata, solution, writeFiles = TRUE)
#' # updating sampling frame with new strata labels
#' data(swissframe)
#' framenew <- updateFrame(frame=swissframe, newstrata=newstrata, writeFiles = TRUE)
#' }
#'
#' @keywords survey
#' @export
updateFrame <- function(frame, 
                        newstrata, 
                        writeFiles = FALSE,
                        outputFolder = file.path(getwd(),"output")) {
  if (writeFiles == TRUE) {
    # if(dir.exists(outputFolder)){
    #   warning("Folder ", outputFolder," already existed and has been deleted.")
    #   unlink(outputFolder)
    # } 
    if(!dir.exists(outputFolder)) dir.create(outputFolder)
  }
    colnames(frame) <- toupper(colnames(frame))
    colnames(newstrata) <- toupper(colnames(newstrata))
    nvarX <- length(grep("X", names(newstrata)))
    if (nvarX == 1) 
        frame$STRATUM <- frame$X1
    if (nvarX > 1) {
        stmt <- NULL
        stmt <- "frame$STRATUM <- paste("
        for (i in 1:(nvarX - 1)) {
            if (i > 0) 
                stmt <- paste(stmt, "frame$X", i, ",", sep = "")
        }
        stmt <- paste(stmt, "frame$X", nvarX, ",sep='*')", sep = "")
        eval(parse(text = stmt))
    }
    newstrata$STRATUM <- as.character(newstrata$STRATUM)
    labels <- as.data.frame(cbind(newstrata$STRATUM, newstrata$DOM1, 
        newstrata$LABEL),stringsAsFactors = TRUE)
    colnames(labels) <- c("STRATUM", "DOMAINVALUE", "LABEL")
    frame$DOMAINVALUE <- as.factor(frame$DOMAINVALUE)  # new line
    framenew <- merge(frame, labels, by = c("DOMAINVALUE", "STRATUM"))
	framenew$STRATUM <- as.factor(framenew$STRATUM)
	framenew$LABEL <- as.character(framenew$LABEL)
	framenew$LABEL <- as.integer(framenew$LABEL)
    colnames(framenew) <- toupper(colnames(framenew))
	if (writeFiles == TRUE) {
		write.table(framenew, file.path(outputFolder,"framenew.txt"), row.names = FALSE, 
			sep = "\t", quote = FALSE)
	}
    # if (writeFiles == TRUE) {
    #   setwd(dire)
    # }
    return(framenew)
}
