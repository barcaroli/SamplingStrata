# ==============================================================================
# File:      selectSampleSytematic.R
# Functions: selectSampleSystematic() [exported]
#
# Purpose:
#   Selects a stratified systematic sample from the frame labelled with the
#   optimized strata: within each stratum the units, optionally sorted by a
#   given variable, are selected starting from a random point with a step equal
#   to the inverse of the sampling rate; returns the sample with weights.
# ==============================================================================

# --------------------------------------------------------------------------
# Function for selecting a sample on the basis of the
# result of optimal stratification and allocation 
# (systematic selection)
# Authors: Giulio Barcaroli 
# Date: 1 April 2018
# --------------------------------------------------------------------------
#' @title Selection of a stratified sample from the frame with systematic method
#'
#' @description
#' Once optimal stratification has been obtained,
#' and a new frame has been
#' built by assigning to the units of the old one the new stratum labels, it is possible to select
#' a stratified sample from the frame
#' with the systematic method, that is a selection that begins selecting the first
#' unit by an initial ramndomly chosen starting point, and proceeding in selecting other units
#' by adding an interval that is the inverse of the sampling rate in the stratum.
#'
#' This selection method can be useful if associated to a particular ordering of the selection
#' frame, where the ordering variable(s) can be considered as additional
#' stratum variable(s).
#'
#' The result of the execution of "selectSampleSystematic" function is a dataframe containing selected
#' units, with the probabilities of inclusion.
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
#'   Name of stratum variable must be 'strato'.
#' @param outstrata This is the (mandatory) dataframe containing the information related to resulting
#'   stratification obtained by the execution of "optimizeStrata" function.
#'   Name of stratum variable must be 'strato'.
#' @param sortvariable This is the name of the variable to be used as ordering variable inside each stratum before proceeding to the systematic selection. It must be previously added to the selection frame. Default is NULL.
#' @param writeFiles Indicates if at the end of the processing the resulting strata will be outputted in a delimited file.
#'   Default is "FALSE".
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
#' newstrata <- updateStrata(swissstrata, solution)
#' # updating sampling frame with new strata labels
#' data(swissframe)
#' framenew <- updateFrame(frame=swissframe,newstrata=newstrata)
#' # adding variable "POPTOT" to framenew
#' data("swissmunicipalities")
#' framenew <- merge(framenew,swissmunicipalities[,c("REG","Nom","POPTOT")],
#'                   by.x=c("REG","ID"),by.y=c("REG","Nom"))
#' # selection of sample with systematic method
#' sample <- selectSampleSystematic(frame=framenew,
#'                 outstrata=solution$aggr_strata,
#'                 sortvariable="POPTOT")
#' head(sample)
#' }
#'
#' @keywords survey
#' @export
selectSampleSystematic <- function(frame, 
                               outstrata, 
                               sortvariable = NULL,
                               writeFiles = FALSE,
                               verbatim = TRUE) {
  strata.sample.systematic <- function(frame, strata, nh, sortvariable, repl) {
    selected <- NULL
    WEIGHTS <- NULL
    # for (i in (1:length(nh))) {
    j <- 0
    for (i in (c(unique(as.numeric(frame$STRATO))))) {
      j <- j+1
      framestrat <- frame[frame[,strata] == i,]
      if (!is.null(sortvariable)) {
        framestrat <- framestrat[order(framestrat[,c(sortvariable)]),]
      }
      # systematic selection with a (possibly fractional) interval:
      # random start uniform in (0, step], then units ceiling(start + k*step),
      # k = 0, ..., n-1. Exactly n units are selected and every unit has
      # inclusion probability n / Nh.
      Nh <- nrow(framestrat)
      n <- nh[j]
      if (n >= Nh) {
        s <- seq_len(Nh)
      } else if (n <= 0) {
        s <- integer(0)
      } else {
        step <- Nh / n
        start <- step - runif(1, 0, step)   # in (0, step]
        s <- ceiling(start + (0:(n - 1)) * step)
      }
      sel <- framestrat$ID[s]
      wgt <- rep(nrow(framestrat)/length(sel),length(sel))
      selected <- c(selected,sel)
      WEIGHTS <- c(WEIGHTS,wgt)
    }
    attr(selected, "WEIGHTS") <- WEIGHTS
    selected
  }
    colnames(frame) <- toupper(colnames(frame))
    if(is.factor(frame$ID)) frame$ID <- as.character(frame$ID)
    colnames(outstrata) <- toupper(colnames(outstrata))
    if (!is.null(sortvariable)) {
      sortvariable <- toupper(sortvariable)
      if (!(sortvariable %in% colnames(frame))) {
        stop("Sort variable ", sortvariable, " not in frame")
      }
    }
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
			s <- strata.sample.systematic(listaord, c("STRATO"), strataord$SOLUZ, 
				sortvariable, repl = FALSE)
			samp <- data.frame(listaord[listaord$ID %in% s, ], WEIGHTS = attr(s, "WEIGHTS"), stringsAsFactors = TRUE)
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
		s <- strata.sample.systematic(listaord, c("STRATO"), strataord$SOLUZ, 
				sortvariable, repl = FALSE)
		samp <- data.frame(listaord[listaord$ID %in% s, ], WEIGHTS = attr(s, "WEIGHTS"), stringsAsFactors = TRUE)
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
        write.table(samptot, "sample.csv", sep = ",", row.names = FALSE, 
            col.names = TRUE, quote = FALSE)
    if (writeFiles == TRUE) 
        write.table(chktot, "sampling_check.csv", sep = ",", 
            row.names = FALSE, col.names = TRUE, quote = FALSE)
    outstrata$FPC <- outstrata$SOLUZ/outstrata$N
	fpc <- outstrata[, c("DOM1","STRATO","FPC")]
	samptot <- merge(samptot, fpc, by.x = c("DOMAINVALUE","STRATO"),by.y=c("DOM1","STRATO"),all.x=TRUE)
    return(samptot)
}
