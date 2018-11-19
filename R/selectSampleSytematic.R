# --------------------------------------------------------------------------
# Function for selecting a sample on the basis of the
# result of optimal stratification and allocation 
# (systematic selection)
# Authors: Giulio Barcaroli 
# Date: 1 April 2018
# --------------------------------------------------------------------------
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
      step <- nrow(framestrat) / nh[j]
      start <- sample(step,1)
      s <- round(seq(start,nrow(framestrat),step),0)
      if (length(s) < nh[j]) {
        s <- c(1,s)
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
        cat("\n Sort variable not in frame")
        stop
      }
    }
    outstrata$SOLUZ <- round(outstrata$SOLUZ)  # rounding of allocation numbers
    numdom <- length(levels(droplevels(as.factor(frame$DOMAINVALUE))))
    samptot <- NULL
    chktot <- NULL
    # begin domains cycle
	if (numdom > 1) {
		for (d in (1:numdom)) {
			domframe <- frame[frame$DOMAINVALUE == d, ]
			domstrata <- outstrata[outstrata$DOM1 == d, ]
			strataord <- domstrata[order(as.numeric(domstrata$STRATO)), ]
			lista <- domframe
			lista$STRATO <- lista$LABEL
			listaord <- lista[order(as.numeric(lista$STRATO)), ]
			s <- strata.sample.systematic(listaord, c("STRATO"), strataord$SOLUZ, 
				sortvariable, repl = FALSE)
			samp <- data.frame(listaord[listaord$ID %in% s, ], WEIGHTS = attr(s, "WEIGHTS"))
			samptot <- rbind(samptot, samp)
			chk <- data.frame(DOMAINVALUE = d, STRATO = strataord$STRATO, 
				Nh_frame = as.vector(table(listaord$STRATO)), Nh_strata = strataord$N, 
				planned_units = strataord$SOLUZ, selected_units = as.vector(table(samp$STRATO)), 
				sum_of_wgts = tapply(samp$WEIGHTS, samp$STRATO, sum))
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
		samp <- data.frame(listaord[listaord$ID %in% s, ], WEIGHTS = attr(s, "WEIGHTS"))
		samptot <- rbind(samptot, samp)
		chk <- data.frame(DOMAINVALUE = strataord$DOM1, STRATO = strataord$STRATO, 
				Nh_frame = as.vector(table(listaord$STRATO)), Nh_strata = strataord$N, 
				planned_units = strataord$SOLUZ, selected_units = as.vector(table(samp$STRATO)), 
				sum_of_wgts = tapply(samp$WEIGHTS, samp$STRATO, sum))
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
