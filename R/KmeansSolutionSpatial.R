KmeansSolutionSpatial <- function (frame, 
                                   fitting = 1, 
                                   range = c(0), 
                                   kappa = 3,
                                   errors, 
                                   nstrata = NA, 
                                   minnumstrat = 2, 
                                   maxclusters = NA, 
                                   showPlot = TRUE) 
{
  cv <- errors
  nvariables <- ncol(cv) - 2
  stmt1 <- "solution <- kmeans(frame[,"
  stmt2 <- "c("
  for (i in (1:nvariables)) {
    if (i < nvariables) 
      stmt2 <- paste(stmt2, "'Y", i, "',", sep = "")
    if (i == nvariables) 
      stmt2 <- paste(stmt2, "'Y", i, "')", sep = "")
  }
  suggestions <- NULL
  domainvalue <- NULL
  id <- NULL
  ndom <- nrow(cv)
  solution <- NULL
  best <- rep(0, ndom)
  best_num_strata <- rep(0, ndom)
  for (k in (1:nrow(cv))) {
    stratacorr <- frame[frame$domainvalue == k, ]
    errorscorr <- errors[errors$domainvalue == k, ]
    aggr <- aggrStrataSpatial(dataset = frame, 
                              fitting = fitting, 
                              range = range, 
                              kappa = kappa, 
                              vett = rep(1, nrow(frame[frame$domainvalue == k,])), 
                              dominio = k)
    v <- bethel(aggr, errorscorr, minnumstrat = minnumstrat)
    sum(v)
    best[k] <- sum(v)
    if (is.na(maxclusters)) {
      times <- round(nrow(stratacorr) * 0.5, 0)
    }
    if (!is.na(maxclusters)) {
      times <- min(maxclusters, round(nrow(stratacorr) * 
                                        0.5, 0))
    }
    if (showPlot == TRUE) {
      plot(1, sum(v), xlim = c(1, times), ylim = c(0, 
                                                   1.5 * sum(v)), type = "p", ylab = "Sample size", 
           xlab = "Number of clusters")
      tit <- paste("title('kmeans clustering in domain ", 
                   k, "')", sep = "")
      eval(parse(text = tit))
    }
    bestsolution <- NULL
    if (is.na(nstrata)) {
      min = 2
      max = times
    }
    if (!is.na(nstrata)) {
      min = nstrata
      max = nstrata
    }
    for (i in min:max) {
      stmt3 <- paste("],", i, ")$cluster", sep = "")
      stmt <- paste(stmt1, stmt2, stmt3, sep = "")
      eval(parse(text = stmt))
      aggr <- aggrStrataSpatial(dataset = stratacorr, 
                                fitting = fitting, 
                                range = range, 
                                kappa = kappa,
                                vett = solution[frame$domainvalue == k], 
                                dominio = k)
      v <- bethel(aggr, errorscorr, minnumstrat = minnumstrat)
      if (showPlot == TRUE) 
        points(i, sum(v))
      if (sum(v) <= best[k]) {
        bestsolution <- solution[frame$domainvalue == k]
        best_num_strata[k] <- i
        best[k] <- sum(v)
      }
    }
    bestsolution <- as.factor(bestsolution)
    levels(bestsolution) <- c(1:length(levels(bestsolution)))
    suggestions <- c(suggestions, bestsolution)
    domainvalue <- c(domainvalue, rep(k, nrow(stratacorr)))
    id <- c(id,as.character(frame$id[frame$domainvalue == k]))
  }
  cat("\n-----------------")
  cat("\n Kmeans solution ")
  cat("\n-----------------")
  for (i in c(1:ndom)) {
    cat("\n *** Domain: ", i, " ***")
    cat("\n Number of strata: ", best_num_strata[i])
    cat("\n Sample size     : ", best[i])
  }
  solutionKmean <- as.data.frame(cbind(id,suggestions,domainvalue),stringsAsFactors = TRUE)
  solutionKmean$domainvalue <- as.integer(solutionKmean$domainvalue)
  return(solutionKmean)
}
