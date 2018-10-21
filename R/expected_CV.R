expected_CV <- function (strata) {
  if (is.null(strata$SOLUZ)) stop("There is no allocation of units in strata")
  ndom <- length(unique(strata$DOM1))
  nvars <- (ncol(strata) - 6) / 2
  cv <- matrix(NA,nrow=ndom,ncol=nvars)
  colnames(cv) <- paste("cv(Y",c(1:nvars),")",sep="")
  rownames(cv) <- paste("DOM",c(1:ndom),sep="")
  for (i in (1:ndom)) {
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
      cv[i,j] <- CV 
    }
  }
  cv <- round(cv,3)
  cv
}