adjustSize <- function (size, strata, cens=NULL, minnumstr = 2) 
{ 
  cens <- strata[strata$CENS == 1,]
  strata <- strata[strata$CENS == 0,]
  if (nrow(cens) > 0) size <- size - sum(cens$N)
  
  diff <- size - sum(strata$SOLUZ)
  diffrel <- diff / size
  
  round(sum(strata$SOLUZ) * 1/(1-diffrel)) + sum(cens$N)
  olddiff = 0
  diff <- round(size - sum(strata$SOLUZ))
  if (diff > 0) {
    while ( diff != olddiff) {
      olddiff = diff
      diffrel <- diff/size
      for (i in (1:nrow(strata))) {
        strata$SOLUZ[i] <- round(strata$SOLUZ[i] * 1/(1-diffrel))
        if (strata$SOLUZ[i] > strata$N[i]) strata$SOLUZ[i] <- strata$N[i]
      }
      diff <- round(size - sum(strata$SOLUZ))
      cat("\n",sum(strata$SOLUZ))
    }
  }
  
  if (diff < 0) {
    while ( diff != olddiff) {
      olddiff = diff
      diffrel <- diff/size
      for (i in (1:nrow(strata))) {
        strata$SOLUZ[i] <- round(strata$SOLUZ[i] * 1/(1-diffrel))
        if (strata$N[i] >= minnumstr & strata$SOLUZ[i] < minnumstr) strata$SOLUZ[i] <- minnumstr
        if (strata$N[i] == 1 & strata$SOLUZ[i] < 1) strata$SOLUZ[i] <- 1
      }
      diff <- round(size - sum(strata$SOLUZ))
      cat("\n",sum(strata$SOLUZ))
    }
  }
  if (nrow(cens) > 0) strata <- rbind(strata,cens)
  cat("\n Final adjusted size: ",sum(strata$SOLUZ))
  return(strata)
}