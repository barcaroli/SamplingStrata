summaryStrata <- function (x, outstrata) 
{
  colnames(x) <- toupper(colnames(x))
  domains <- length(levels(as.factor(x$DOMAINVALUE)))
  nvars <- length(grep("M", colnames(outstrata)))-1
  vars <- colnames(x)[grep("X", colnames(x))]
  outstrata <- outstrata[order(paste(as.numeric(outstrata$DOM1), 
                                     as.numeric(outstrata$STRATO))), ]
  out <- NULL
  out$Domain <- outstrata$DOM1
  out$Stratum <- outstrata$STRATO
  out$Population <- outstrata$N
  out$Allocation <- round(outstrata$SOLUZ)
  out$"Sampling rate" <- outstrata$SOLUZ/outstrata$N
  out <- as.data.frame(out)
  colnames(out) <- c("Domain", "Stratum", "Population", "Allocation", 
                     "SamplingRate")
  for (i in 1:nvars) {
    stringa <- paste("out$bounds_X", i, " <- rep(NA,nrow(out))", 
                     sep = "")
    eval(parse(text = stringa))
    lab <- paste("Bounds_X", i, "", sep = "")
    eval(parse(text = stringa))
    colnames(out)[ncol(out)] <- lab
  }
  for (j in 1:domains) {
    xdom <- x[x$DOMAINVALUE == j, ]
    nstrata <- length(unique(xdom$LABEL))
    for (i in 1:nvars) {
      stringa <- paste("x", i, "_min <- tapply(xdom$", 
                       vars[i], ",xdom$LABEL,min)", sep = "")
      eval(parse(text = stringa))
      stringa <- paste("x", i, "_max <- tapply(xdom$", 
                       vars[i], ",xdom$LABEL,max)", sep = "")
      eval(parse(text = stringa))
      stringa <- paste("xcuts <- c(c(x", i, "_min[2:(length(x", 
                       i, "_min)-1)]),x", i, "_max[length(x", i, "_min)],x", 
                       i, "_max[length(x", i, "_max)])", sep = "")
      eval(parse(text = stringa))
      stringa <- paste("out$Bounds_X", i, "[out$Domain == j] <- paste(x", 
                       i, "_min,x", i, "_max,sep=' - ')", sep = "")
      eval(parse(text = stringa))
    }
  }
  return(out)
}
