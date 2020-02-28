# Function to prepare suggestions from KmeansSolution2
# or from KmeansSolutionSpatial for optimizeStrata2
# or optimizeStrataSpatial

prepareSuggestion <- function(kmean=kmean,
                              frame=frame,
                              nstrat=nstrat) {
  # kmean$id <- row.names(kmean)
  frame1 <- frame
  nvarX <- sum(grepl("X",colnames(frame1)))
  for (i in (1:nvarX)) {
    st <- paste("frame1$X",i," <- NULL",sep="")
    eval(parse(text=st))
  }
  frame1 <- merge(frame1,kmean[,c("id","suggestions")],by=c("id"))
  frame1$X1 <- frame1$suggestions
  strataKm <- buildStrataDF(frame1,progress=FALSE)
  # strataKm$SOLUZ <- bethel(strataKm,cv[1,])
  strataKm$SOLUZ <- 1
  framenew <- frame
  framenew <- merge(framenew,kmean[,c("id","suggestions")],by=c("id"))
  framenew$LABEL <- frame1$suggestions
  ss <- summaryStrata(framenew,strataKm,progress=FALSE)
  ndom <- length(unique(ss$Domain))
  nvarX <- length(grep("X",colnames(ss)))/2
  ss <- ss[order(ss$Lower_X1),]
  ss <- ss[order(ss$Domain),]
  nr <- 0
  for (k in c(1:ndom)) {
    nr <- nr + (nstrat[k]-1) * nvarX
  }
  sugg <- NULL
  sugg$domainvalue <- rep(NA,nr)
  sugg$suggestions1 <- rep(NA,nr)
  sugg$suggestions2 <- rep(NA,nr)
  sugg <- as.data.frame(sugg,stringsAsFactors = TRUE)
  for (d in (1:length(unique(ss$Domain)))) {
    jmin <- min(which(is.na(sugg$suggestions1)))
    jmax <- jmin - 1 + (nstrat[d] - 1) * nvarX
    sugg$domainvalue[jmin:jmax] <- d
    for (i in (1:nvarX)) {
      kmin <- jmin + (i-1)*(nstrat[d]-1)
      kmax <- jmin - 1 + i*(nstrat[d] - 1) 
      stmt <- paste("sugg$suggestions1[kmin:kmax] <- (ss$Lower_X",i,"[order(ss$Lower_X",i,")]/max(frame$X",i,"))[2:(nstrat[d])]",sep="")
      eval(parse(text=stmt))
      stmt <- paste("sugg$suggestions2[kmin:kmax] <- (ss$Upper_X",i,"[order(ss$Upper_X",i,")]/max(frame$X",i,"))[1:(nstrat[d]-1)]",sep="")
      eval(parse(text=stmt))
    }
  }
  return(sugg)
}
