summaryStrata <- function (x, 
                           outstrata, 
                           progress = TRUE, 
                           writeFiles = FALSE) 
{
  colnames(x) <- toupper(colnames(x))
  domains <- length(levels(as.factor(x$DOMAINVALUE)))
  doms <- as.numeric(levels(as.factor(x$DOMAINVALUE)))
  nvars <- length(grep("X", colnames(x))) 
  vars <- colnames(x)[grep("X", colnames(x))]
  outstrata <- outstrata[order(as.numeric(outstrata$DOM1), 
                               as.numeric(outstrata$STRATO)), ]
  out <- NULL
  out$Domain <- outstrata$DOM1
  out$Stratum <- outstrata$STRATO
  out$Population <- outstrata$N
  out$Allocation <- round(outstrata$SOLUZ)
  out$"Sampling rate" <- round(outstrata$SOLUZ/outstrata$N, 
                               6)
  out <- as.data.frame(out,stringsAsFactors = TRUE)
  colnames(out) <- c("Domain", "Stratum", "Population", "Allocation", 
                     "SamplingRate")
  for (i in 1:nvars) {
    stringa <- paste("out$min_X", i, " <- rep(NA,nrow(out))", 
                     sep = "")
    eval(parse(text = stringa))
    stringa <- paste("out$max_X", i, " <- rep(NA,nrow(out))", 
                     sep = "")
    eval(parse(text = stringa))
    lab <- paste("Lower_X", i, "", sep = "")
    eval(parse(text = stringa))
    colnames(out)[ncol(out) - 1] <- lab
    lab <- paste("Upper_X", i, "", sep = "")
    eval(parse(text = stringa))
    colnames(out)[ncol(out)] <- lab
  }
  if (progress == TRUE) 
    pb <- txtProgressBar(min = 0, max = domains, style = 3)
  for (j in doms) {
    if (progress == TRUE) 
      Sys.sleep(0.1)
    if (progress == TRUE) 
      setTxtProgressBar(pb, j)
    xdom <- x[x$DOMAINVALUE == j, ]
    nstrata <- length(unique(xdom$LABEL))
    for (i in 1:nvars) {
      stringa <- paste("x", i, "_min <- tapply(xdom$", 
                       vars[i], ",xdom$LABEL,min)", sep = "")
      eval(parse(text = stringa))
      stringa <- paste("x", i, "_max <- tapply(xdom$", 
                       vars[i], ",xdom$LABEL,max)", sep = "")
      eval(parse(text = stringa))
      stringa <- paste("out$Lower_X", i, "[out$Domain == j] <- x", 
                       i, "_min", sep = "")
      eval(parse(text = stringa))
      stringa <- paste("out$Upper_X", i, "[out$Domain == j] <- x", 
                       i, "_max", sep = "")
      eval(parse(text = stringa))
    }
  }
  if (progress == TRUE) 
    close(pb)
  if (writeFiles == TRUE) 
    write.table(out, "strata_structure.csv", sep = ";", 
                row.names = FALSE)
  return(out)
}
