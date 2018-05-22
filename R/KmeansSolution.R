KmeansSolution <- function(strata,
                             errors,
                             nstrata = NA,
                             maxclusters = NA,
                             showPlot = FALSE) {
  nvariables <- ncol(errors) - 2
  stmt1 <- "solution <- kmeans(stratacorr[,"
  stmt2 <- "c("
  for (i in (1:nvariables)) {
    if (i < nvariables) stmt2 <- paste(stmt2,"'M",i,"',",sep="")
    if (i == nvariables) stmt2 <- paste(stmt2,"'M",i,"')",sep="")
  }
  suggestions <- NULL
  domainvalue <- NULL
  solution <- NULL
    for (k in (unique(strata$DOM1))) {
      stratacorr <- strata[strata$DOM1 == k,]
      errorscorr <- errors[errors$domainvalue == k,]
      aggr <- aggrStrata(strata=stratacorr,
                         nvar=nvariables,
                         censiti=0,
                         vett=rep(1,nrow(stratacorr)),
                         dominio=k)
      v <- bethel(aggr, errorscorr)
      sum(v)
      best <- sum(v)
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
        v <- bethel(aggr, errorscorr, minnumstrat=2)
        # cat("\n",sum(v))
        if (showPlot == TRUE) points(i,sum(v))
        if (sum(v) <= best) {
          bestsolution <- solution
          best_num_strata <- i
          best <- sum(v)
        }
      }
      suggestions <- c(suggestions,bestsolution)
      domainvalue <- c(domainvalue,rep(k,nrow(stratacorr)))
  }
  solutionKmean <- as.data.frame(cbind(suggestions,domainvalue))
  cat("\n Kmeans solution")
  cat("\n Number of strata: ",best_num_strata)
  cat("\n Sample size     : ",best)
  solutionKmean
}