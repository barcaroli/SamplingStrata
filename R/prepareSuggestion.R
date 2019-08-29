# Function to prepare suggestions from KmeansSolution2
# or from KmeansSolutionSpatial for optimizeStrata2
# or optimizeStrataSpatial

prepareSuggestion <- function(ss) {
  ndom <- length(unique(ss$Domain))
  nvarX <- length(grep("X",colnames(ss)))/2
  nstrat <- nrow(ss[ss$Domain == 1,])
  ss <- ss[order(ss$Lower_X1),]
  sugg <- NULL
  sugg$domainvalue <- rep(NA,ndom*(nstrat-1)*nvarX)
  sugg$suggestions1 <- rep(NA,ndom*(nstrat-1)*nvarX)
  sugg$suggestions2 <- rep(NA,ndom*(nstrat-1)*nvarX)
  sugg <- as.data.frame(sugg)
  for (d in (1:length(unique(ss$Domain)))) {
    jmax <- d*(nstrat-1)*nvarX
    jmin <- jmax - (nstrat-1)*nvarX + 1
    sugg$domainvalue[jmin:jmax] <- d
    for (i in (1:nvarX)) {
      kmin <- jmin + (i-1)*(nstrat-1)
      kmax <- jmin + i*(nstrat-1) - 1
      stmt <- paste("sugg$suggestions1[kmin:kmax] <- (ss$Lower_X",i,"[order(ss$Lower_X",i,")]/max(frame$X",i,"))[2:(nstrat)]",sep="")
      eval(parse(text=stmt))
      stmt <- paste("sugg$suggestions2[kmin:kmax] <- (ss$Upper_X",i,"[order(ss$Upper_X",i,")]/max(frame$X",i,"))[1:(nstrat-1)]",sep="")
      eval(parse(text=stmt))
    }
  }
  return(sugg)
}
