<<<<<<< HEAD
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
=======
prepareSuggestion <- function(frame,model,solution_indices) {
  colnames(frame) <- toupper(colnames(frame))
  framereduced <- frame
  framereduced[,c(grep("X",colnames(framereduced)))] <- NULL
  framereduced$X1 <- solution_indices
  strata <- buildStrataDF2(dataset=framereduced,model=model,progress=FALSE,verbose=FALSE)
  strata$SOLUZ <- rep(2,nrow(strata))
  frame$LABEL <- solution_indices
  ss <- summaryStrata(frame,strata)
  nvarX <- length(grep("X",colnames(ss)))/2
  for (i in (1:nvarX)) {
    stmt <- paste("cuts_X",i," <- NULL",sep="")
    eval(parse(text=stmt))
  }
  suggest <- NULL
  sugg <- NULL
  for (d in (unique(ss$Domain))) {
    sugg$domainvalue <- d
    ssreduced <- ss[ss$Domain==d,]
    nstrat <- nrow(ssreduced)
    stringa <- "sugg$lower <- c("
    for (i in (1:nvarX)) {
      stmt <- paste("cuts_X",i," <- (ssreduced$Lower_X",i,"[order(ssreduced$Lower_X",i,")]/max(ssreduced$Upper_X",i,"))[2:nstrat]",sep="")
      eval(parse(text=stmt))
      stringa <- paste(stringa,"cuts_X",i,sep="")
      if (i < nvarX) stringa <- paste(stringa,",",sep="")
      if (i == nvarX) stringa <- paste(stringa,")",sep="")
    }
    eval(parse(text=stringa))
    stringa <- "sugg$upper <- c("
    for (i in (1:nvarX)) {
      stmt <- paste("cuts_X",i," <- (ssreduced$Upper_X",i,"[order(ssreduced$Upper_X",i,")]/max(ssreduced$Upper_X",i,"))[1:(nstrat-1)]",sep="")
      eval(parse(text=stmt))
      stringa <- paste(stringa,"cuts_X",i,sep="")
      if (i < nvarX) stringa <- paste(stringa,",",sep="")
      if (i == nvarX) stringa <- paste(stringa,")",sep="")
    }
    eval(parse(text=stringa))
    sugg <- as.data.frame(sugg)
    suggest <- rbind(suggest,sugg)
  }
  
  return(suggest)
}
>>>>>>> b2d04e628d4fbd96c8264c709f5ad93d4ad91516
