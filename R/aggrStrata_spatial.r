# ----------------------------------------------------
# Function to aggregate initial atomic strata
# on the basis of the current solution 
# expressed by a vector of integer values
# Version valid for the spatial case
# Author: Giulio Barcaroli
# Date: 6 April 2019
# ----------------------------------------------------
aggrStrataSpatial <- function(dataset,
                              fitting=c(1),
                              range=c(0),
                              kappa=3, 
                              vett, 
                              dominio) {

  #---------------------------------------------------------    
  # standard deviation calculated with distances
  stdev <- function(zz, dist, var, STRATO, dataset, fitting, range, kappa) {
    ind <- which(dataset$STRATO == STRATO)
    # differences
    z_z<-zz[ind,ind]
    # distances
    dist <- dist[ind,ind]
    # variances
    var <- var[ind]
    
    if (length(ind) > 1) {
      somma_coppie_var <- as.matrix(outer(var,var,"+"))
      spatial_correlation <- (1 - (exp(-kappa*dist/range)))
    }
    if (length(ind) <= 1) {
      somma_coppie_var <- 0
      spatial_correlation <- 0
    }
    # variance in the stratum
    D2 <- z_z/fitting + somma_coppie_var * spatial_correlation
    var_strato <- sum(D2) / (2*length(ind)^2)
    # standard deviation
    if (var_strato < 0) var_strato <- 0
    sd_strato <- sqrt(var_strato)
    return(sd_strato)
  }
  #--------------------------------------------------------- 
  colnames(dataset) <- toupper(colnames(dataset))
  dataset <- dataset[dataset$DOMAINVALUE == dominio,]
  dist <- sqrt((outer(dataset$LON,dataset$LON,"-"))^2+(outer(dataset$LAT,dataset$LAT,"-"))^2)
  nvarY <- length(grep("Y", names(dataset)))
  stratatot <- NULL
  dataset$STRATO <- as.factor(vett)
  for (j in as.numeric(levels(dataset$STRATO))) {
    strato <- NULL
    strato$stratum <- j
    strato$N <- nrow(dataset[dataset$STRATO == j,])
    strato$COST <- 1
    strato$CENS <- 0
    strato$DOM1 <- dominio
    strato <- as.data.frame(strato, stringsAsFactors = TRUE)
    rng <- NULL
    fit <- NULL
    zz <- NULL
    for (i in 1:nvarY) {
      stmt <- paste("var <- dataset$VAR",i,sep="")
      eval(parse(text = stmt))
      stmt <- paste("strato$M", i, " <- mean(dataset$Y",i,"[dataset$STRATO == ",j,"])",sep = "")
      eval(parse(text = stmt))
      stmt <- paste("zz <- outer(dataset$Y",i,",dataset$Y",i,",'-')^2",sep="")
      eval(parse(text = stmt))
      stmt <- paste("rng <- range[",i,"]",sep="")
      eval(parse(text = stmt))
      stmt <- paste("fit <- fitting[",i,"]",sep="")
      eval(parse(text = stmt))
      sd <- stdev(zz,dist,var,j,dataset,fit,rng,kappa)
      stmt <- paste("strato$S",i," <- sd",sep="")
      eval(parse(text=stmt))
    }
    stratatot <- rbind(stratatot, strato)
    }
  return(stratatot)
}