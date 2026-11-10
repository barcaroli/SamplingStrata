# ==============================================================================
# File:      selectSample.R
# Functions: selectSample() [exported]
#
# Purpose:
#   Selects a stratified simple random sample without replacement (srswor) from
#   the frame labelled with the optimized strata, according to the allocation
#   in the aggregated strata; returns the sample with design weights and,
#   optionally, writes sample and frame to files.
# ==============================================================================

#
# --------------------------------------------------------------------------
# Function for selecting a sample on the basis of the
# result of optimal stratification and allocation 
# Authors: Giulio Barcaroli 
# with a contribution from Diego Zardetto
# Date: 4 January 2012
# Last update: 12 January 2023
# --------------------------------------------------------------------------
#' @title Selection of a stratified sample from the frame with srswor method
#'
#' @description
#' Once optimal stratification has been obtained,
#' and a new frame has been
#' built by assigning to the units of the old one the new strata labels, it is possible to select
#' a stratified sample from the frame
#' with the simple random sampling without replacement (srswor) method.
#' The result of the execution of "selectSample" function is a dataframe containing the selected
#' units, with their weights (inverse of the probabilities of inclusion).
#' It is possible to output this dataframe in a .csv file.
#' One more .csv file is produced ("sampling_check"), containing coeherence checks between
#' (a) population in frame strata
#' (b) population in optimised strata
#' (c) planned units to be selected in optimised strata
#' (d) actually selected units
#' (e) sum of weights in each stratum
#'
#' @param frame This is the (mandatory) dataframe containing the sampling frame, as it has been modified
#'   by the execution of the "updateFrame" function.
#' @param outstrata This is the (mandatory) dataframe containing the information related to resulting   stratification obtained by the execution of "optimizeStrata" function. It should coincide with 'solution$aggr_strata'.
#' @param writeFiles Indicates if at the end of the processing the resulting strata will be outputted in a delimited file.
#'   Default is "FALSE".
#' @param outputFolder A path to a directory ot be used for the outputs. Default is "output"
#' @param verbatim Indicates if information on the drawn sample must be printed or not.
#'   Default is "TRUE".
#'
#' @return
#' A dataframe containing the sample
#'
#' @author
#' Giulio Barcaroli with contribution from Diego Zardetto
#'
#' @examples
#' \dontrun{
#' library(SamplingStrata)
#' data(swisserrors)
#' data(swissstrata)
#' # optimisation of sampling strata
#' solution <- optimizeStrata (
#'     errors = swisserrors,
#'     strata = swissstrata
#' )
#' # updating sampling strata with new strata labels
#' newstrata <- updateStrata(swissstrata, solution)
#' # updating sampling frame with new strata labels
#' data(swissframe)
#' framenew <- updateFrame(frame=swissframe,newstrata=newstrata)
#' # selection of sample
#' sample <- selectSample(frame=framenew,outstrata=solution$aggr_strata)
#' head(sample)
#' }
#'
#' @keywords survey
#' @export
selectSample <- function(frame, 
                         outstrata, 
                         writeFiles = FALSE,
                         verbatim=TRUE,
                         outputFolder = file.path(getwd(),"output")) {
  if (writeFiles == TRUE) {
    # if(dir.exists(outputFolder)){
    #   warning("Folder ", outputFolder," already existed and has been deleted.")
    #   unlink(outputFolder)
    # } 
    if(!dir.exists(outputFolder)) dir.create(outputFolder)
  }
    strata.sample <- function(frame, strata, nh, repl) {
        stratodist <- table(frame[, strata])
        stratocum <- c(0, cumsum(stratodist))
        s.stratocum <- c(0, cumsum(nh))
        permuta <- rep(NA, sum(nh))
        sapply(1:length(nh), function(i) {
            permuta[(s.stratocum[i] + 1):s.stratocum[i + 1]] <<- (stratocum[i] + 
                sample(stratodist[i], nh[i], replace = repl))
        })
        WEIGHTS <- rep(if (repl == FALSE) stratodist/nh else 1/(1 - 
            (1 - 1/stratodist)^nh), nh)
        attr(permuta, "WEIGHTS") <- WEIGHTS
        permuta
    }
    colnames(frame) <- toupper(colnames(frame))
    colnames(outstrata) <- toupper(colnames(outstrata))
    outstrata$SOLUZ <- round(outstrata$SOLUZ)  # rounding of allocation numbers
    # numdom <- length(levels(droplevels(as.factor(frame$DOMAINVALUE))))
    frame$DOMAINVALUE <- factor(frame$DOMAINVALUE)
    numdom <- length(levels(frame$DOMAINVALUE))
    samptot <- NULL
    chktot <- NULL
    # begin domains cycle
	if (numdom > 1) {
	  # for (d in (1:numdom) {
		for (d in (levels(frame$DOMAINVALUE))) {
			domframe <- frame[frame$DOMAINVALUE == d, ]
			domstrata <- outstrata[outstrata$DOM1 == d, ]
			strataord <- domstrata[order(as.numeric(domstrata$STRATO)), ]
			lista <- domframe
			lista$STRATO <- lista$LABEL
			listaord <- lista[order(as.numeric(lista$STRATO)), ]
			s <- strata.sample(listaord, c("STRATO"), strataord$SOLUZ, 
				repl = FALSE)
			samp <- data.frame(listaord[s, ], WEIGHTS = attr(s, "WEIGHTS"), stringsAsFactors = TRUE)
			samptot <- rbind(samptot, samp)
			chk <- data.frame(DOMAINVALUE = d, STRATO = strataord$STRATO, 
				Nh_frame = as.vector(table(listaord$STRATO)), Nh_strata = strataord$N, 
				planned_units = strataord$SOLUZ, selected_units = as.vector(table(samp$STRATO)), 
				sum_of_wgts = tapply(samp$WEIGHTS, samp$STRATO, sum), stringsAsFactors = TRUE)
			chktot <- rbind(chktot, chk)
		}  # end domain cycle
	}
	if (numdom == 1) {
		domframe <- frame
		domstrata <- outstrata
		strataord <- domstrata[order(as.numeric(domstrata$STRATO)), ]
		lista <- domframe
		lista$STRATO <- lista$LABEL
		listaord <- lista[order(lista$STRATO), ]
		s <- strata.sample(listaord, c("STRATO"), strataord$SOLUZ, 
				repl = FALSE)
		samp <- data.frame(listaord[s, ], WEIGHTS = attr(s, "WEIGHTS"), stringsAsFactors = TRUE)
		samptot <- rbind(samptot, samp)
		chk <- data.frame(DOMAINVALUE = strataord$DOM1, STRATO = strataord$STRATO, 
				Nh_frame = as.vector(table(listaord$STRATO)), Nh_strata = strataord$N, 
				planned_units = strataord$SOLUZ, selected_units = as.vector(table(samp$STRATO)), 
				sum_of_wgts = tapply(samp$WEIGHTS, samp$STRATO, sum), stringsAsFactors = TRUE)
		chktot <- rbind(chktot, chk)
	}
    colnames(samptot) <- toupper(colnames(samptot))
    colnames(chktot) <- toupper(colnames(chktot))
    cens <- sum((chktot$NH_STRATA == chktot$PLANNED_UNITS) == 
        TRUE)
    cens.units <- sum(chktot$PLANNED_UNITS[chktot$NH_STRATA == 
        chktot$PLANNED_UNITS])
	if (verbatim == TRUE) { 
		cat("\n*** Sample has been drawn successfully ***")
		cat("\n", nrow(samptot), " units have been selected from ", 
			nrow(outstrata), " strata\n")
		if (cens > 0) {
			cat("\n==> There have been ", cens, " take-all strata ")
			cat("\nfrom which have been selected ", cens.units, "units\n")
		}
	}
    if (writeFiles == TRUE) 
        write.table(samptot, file.path(outputFolder,"sample.csv"), sep = ",", row.names = FALSE, 
            col.names = TRUE, quote = FALSE)
    if (writeFiles == TRUE) 
        write.table(chktot, file.path(outputFolder,"sampling_check.csv"), sep = ",", 
            row.names = FALSE, col.names = TRUE, quote = FALSE)
    outstrata$FPC <- outstrata$SOLUZ/outstrata$N
	fpc <- outstrata[, c("DOM1","STRATO","FPC")]
	samptot <- merge(samptot, fpc, by.x = c("DOMAINVALUE","STRATO"),by.y=c("DOM1","STRATO"),all.x=TRUE)
    return(samptot)
}
