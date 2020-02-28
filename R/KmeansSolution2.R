KmeansSolution2 <- function(frame,
                            model=NULL,
                            errors,
                            nstrata = NA,
                            minnumstrat =2,
                            maxclusters = NA,
                            showPlot = TRUE) {
  colnames(frame) <- toupper(colnames(frame))
  colnames(errors) <- toupper(colnames(errors))
  nvariables <- ncol(errors) - 2
  stmt1 <- "solution <- kmeans(framecorr[,"
  stmt2 <- "c("
  for (i in (1:nvariables)) {
    if (i < nvariables) stmt2 <- paste(stmt2,"'Y",i,"',",sep="")
    if (i == nvariables) stmt2 <- paste(stmt2,"'Y",i,"')",sep="")
  }
  suggestions <- NULL
  domainvalue <- NULL
  id <- NULL
  ndom <- length(unique(frame$DOMAINVALUE))
  solution <- NULL
  best <- NULL
  best_num_strata <- NULL
  for (k in (unique(frame$DOMAINVALUE))) {
    framecorr <- frame[frame$DOMAINVALUE == k,]
    framecorr[,grep("X",colnames(framecorr))] <- NULL
    errorscorr <- errors[errors$DOMAINVALUE == k,]
    stmt3 <- paste("],2)$cluster",sep="")
    stmt <- paste(stmt1,stmt2,stmt3,sep="")
    eval(parse(text=stmt))
    framecorr$X1 <- solution
    aggr <- buildStrataDF(dataset=framecorr,
                           model=NULL,
                           progress=FALSE,
                           verbose=FALSE)
    v <- bethel(aggr, errorscorr, minnumstrat = minnumstrat)

    best[k] <- sum(v)
    if (is.na(maxclusters)) {
      times <- round(nrow(framecorr)*0.5,0)
    }
    if (!is.na(maxclusters)) {
      times <- min(maxclusters,round(nrow(framecorr)*0.5,0))
    }
    if (showPlot == TRUE) {
  	  plot(2,sum(v),xlim=c(1,times),ylim=c(0,1.5*sum(v)),type="p",
         ylab="Sample size",xlab="Number of clusters")
  	  tit <- paste("title('kmeans clustering in domain ",k,"')",sep="")
  	  eval(parse(text=tit))
    }
    bestsolution <- NULL
    if (is.na(nstrata)) {min = 3; max = times }
    if (!is.na(nstrata)) {min = nstrata; max = nstrata }
    for (i in min:max) {
      stmt3 <- paste("],",i,")$cluster",sep="")
      stmt <- paste(stmt1,stmt2,stmt3,sep="")
      eval(parse(text=stmt))
      framecorr$X1 <- solution
      aggr <- buildStrataDF(dataset=framecorr,
                             model=NULL,
                             progress=FALSE,
                             verbose=FALSE)
      v <- bethel(aggr, errorscorr, minnumstrat = minnumstrat)
      # cat("\n",sum(v))
      if (showPlot == TRUE) points(i,sum(v))
      if (sum(v) <= best[k]) {
        bestsolution <- solution
        best_num_strata[k] <- i
        best[k] <- sum(v)
        cat("\nStrata ",i)
        cat("\nSample size ",sum(v))
      }
    }
    bestsolution <- as.factor(bestsolution)
    levels(bestsolution) <- c(1:length(levels(bestsolution)))
    suggestions <- c(suggestions,bestsolution)
    domainvalue <- c(domainvalue,rep(k,nrow(framecorr)))
    id <- c(id,as.character(frame$ID[frame$DOMAINVALUE == k]))
  }
  cat("\n-----------------")
  cat("\n Kmeans solution ")
  cat("\n-----------------")
  for (i in c(1:ndom)) {
    cat("\n *** Domain: ",i," ***")
    cat("\n Number of strata: ",best_num_strata[i])
    cat("\n Sample size     : ",best[i])
  }
  solutionKmean <- as.data.frame(cbind(id,suggestions,domainvalue),stringsAsFactors = TRUE)
  solutionKmean$domainvalue <- as.integer(solutionKmean$domainvalue)
  return(solutionKmean)
}
