KmeansSolution <- function(strata,
                             errors,
                             nstrata = NA,
                             minnumstrat =2,
                             maxclusters = NA,
                             showPlot = TRUE) {
  nvariables <- ncol(errors) - 2
  stmt1 <- "solution <- kmeans(stratacorr[,"
  stmt2 <- "c("
  for (i in (1:nvariables)) {
    if (i < nvariables) stmt2 <- paste(stmt2,"'M",i,"',",sep="")
    if (i == nvariables) stmt2 <- paste(stmt2,"'M",i,"')",sep="")
  }
  suggestions <- NULL
  domainvalue <- NULL
  ndom <- length(unique(strata$DOM1))
  solution <- NULL
  best <- rep(0,ndom)
  best_num_strata <- rep(0,ndom)
    for (k in (unique(strata$DOM1))) {
      stratacorr <- strata[strata$DOM1 == k,]
      errorscorr <- errors[errors$domainvalue == k,]
      aggr <- aggrStrata(strata=stratacorr,
                         nvar=nvariables,
                         censiti=0,
                         vett=rep(1,nrow(stratacorr)),
                         dominio=k)
      v <- bethel(aggr, errorscorr, minnumstrat = minnumstrat)
      sum(v)
      best[k] <- sum(v)
      if (is.na(maxclusters)) {
        times <- round(nrow(stratacorr)*0.5,0)
      }
      if (!is.na(maxclusters)) {
        times <- min(maxclusters,round(nrow(stratacorr)*0.5,0))
      }
      if (showPlot == TRUE) {
    	  plot(1,sum(v),xlim=c(1,times),ylim=c(0,1.5*sum(v)),type="p",
           ylab="Sample size",xlab="Number of clusters")
    	  tit <- paste("title('kmeans clustering in domain ",k,"')",sep="")
    	  eval(parse(text=tit))
      }
      bestsolution <- NULL
      if (is.na(nstrata)) {min = 2; max = times }
      if (!is.na(nstrata)) {min = nstrata; max = nstrata }
      for (i in min:max) {
        stmt3 <- paste("],",i,")$cluster",sep="")
        stmt <- paste(stmt1,stmt2,stmt3,sep="")
        eval(parse(text=stmt))
        aggr <- aggrStrata(strata=stratacorr,
                           nvar=nvariables,
                           censiti=0,
                           vett=solution,
                           dominio=k)
        v <- bethel(aggr, errorscorr, minnumstrat = minnumstrat)
        # cat("\n",sum(v))
        if (showPlot == TRUE) points(i,sum(v))
        if (sum(v) <= best[k]) {
          bestsolution <- solution
          best_num_strata[k] <- i
          best[k] <- sum(v)
        }
      }
      suggestions <- c(suggestions,bestsolution)
      domainvalue <- c(domainvalue,rep(k,nrow(stratacorr)))
  }
  solutionKmean <- as.data.frame(cbind(suggestions,domainvalue),stringsAsFactors = TRUE)
  solutionKmean$domainvalue <- as.integer(solutionKmean$domainvalue)
  totsize <- 0
  cat("\n-------------------")
  cat("\n  Kmeans solution ")
  cat("\n-------------------")
  for (i in c(1:ndom)) {
    cat("\n *** Domain: ", i, " ***")
    cat("\n Number of strata: ", best_num_strata[i])
    cat("\n Sample size     : ", best[i])
    totsize <- totsize + best[i]
  }
  cat("\n-------------------")
  cat("\n Total size: ",totsize)
  cat("\n-------------------")
  solutionKmean
}