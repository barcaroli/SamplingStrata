expected_CV <- function (strata) {
  M_h <- S_h <- NULL
  if (is.null(strata$SOLUZ)) stop("There is no allocation of units in strata")
  ndom <- length(unique(strata$DOM1))
  nvarX <- length(grep("X",colnames(strata)))
  nstr <- length(grep("STR",colnames(strata)))
  nvars <- (ncol(strata) - nstr - 5 - nvarX) / 2
  cv <- matrix(NA,nrow=ndom,ncol=nvars)
  colnames(cv) <- paste("cv(Y",c(1:nvars),")",sep="")
  rownames(cv) <- paste("DOM",c(1:ndom),sep="")
  k<-0
  for (i in (as.numeric(levels(as.factor(strata$DOM1))))) {
    k<-k+1
    stratadom <- strata[strata$DOM1 == i,]
    for (j in 1:nvars) {
      n_h <- stratadom$SOLUZ
      N_h <- stratadom$N
      stmt <- paste("S_h <- stratadom$S",j,sep="")
      eval(parse(text=stmt))
      stmt <- paste("M_h <- stratadom$M",j,sep="")
      eval(parse(text=stmt))
      Y_h <- N_h * M_h
      Var_h <- (N_h^2) * (1 - n_h/N_h) * ((S_h^2)/n_h)
      CV <- sqrt(sum(Var_h)) / sum(Y_h)
      cv[k,j] <- CV 
    }
  }
  cv <- round(cv,7)
  cv
}