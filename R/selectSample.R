#
# --------------------------------------------------------------------------
# Function for selecting a sample on the basis of the
# result of optimal stratification and allocation 
# Authors: Giulio Barcaroli 
# with a contribution from Diego Zardetto
# Date: 4 January 2012
# Last update: 12 January 2023
# --------------------------------------------------------------------------
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
