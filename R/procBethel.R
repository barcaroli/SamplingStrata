#-------------------------------------------------------
# Function to optimize allocation and to select a sample
# with a direct Bethel multivariate allocation
# (without application of genetic algorithm)
#-------------------------------------------------------
procBethel <- function(framesamp,
                        framecens,
                        errors,
                        sampling_method=c("srs","systematic","spatial"),
                        minnumstrat=2) {
  # Sampling frame -----------------------------------------------------
  colnames(framesamp) <- toupper(colnames(framesamp))
  nvarX <- length(grep("X",colnames(framesamp)))
  nvarY <- length(grep("Y",colnames(framesamp)))
  st1 <- "framesamp$STRATUM <- paste(framesamp$DOMAINVALUE,framesamp$X"
  for (i in c(1:nvarX)) {
    if (i < nvarX) st1 <- paste0(st1,i,",framesamp$X")
    if (i == nvarX) st1 <- paste0(st1,i,",sep='*')")
  }
  eval(parse(text=st1))
  framesamp$STRATO <- as.numeric(as.factor(framesamp$STRATUM))
  framesamp$LABEL <- framesamp$STRATO
  # Strata         -----------------------------------------------------
  strata <- buildStrataDF(framesamp,progress=F)
  colnames(strata) <- toupper(colnames(strata))
  st2 <- "strata$STRATUM <- paste(strata$DOM1,strata$X"
  for (i in c(1:nvarX)) {
    if (i < nvarX) st2 <- paste0(st2,i,",strata$X")
    if (i == nvarX) st2 <- paste0(st2,i,",sep='*')")
  }
  eval(parse(text=st2))
  strata$STRATO <- as.numeric(as.factor(strata$STRATUM))
  # Eliminate rare cases
  # for (i in (1:nvarY)) {
  #   st <- paste0("strata$S",i," <- ifelse(strata$M",i," < 0.01,0,strata$S",i,")")
  #   eval(parse(text=st))
  # }

  # Take-all strata ----------------------------------------------------
  if (!is.null(framecens)) {
    cens <- buildStrataDF(framecens,progress=F)
    cens$CENS <- 1
    cens$SOLUZ <- cens$N
    colnames(framecens) <- toupper(colnames(framecens))
    cens$STRATUM <- as.numeric(as.factor(cens$STRATO))
    st1 <- "framecens$STRATUM <- paste(framecens$DOMAINVALUE,framecens$X"
    for (i in c(1:nvarX)) {
      if (i < nvarX) st1 <- paste0(st1,i,",framecens$X")
      if (i == nvarX) st1 <- paste0(st1,i,",sep='*')")
    }
    eval(parse(text=st1))
    framecens$LABEL <- as.numeric(as.factor(framecens$STRATUM))
    framecens$STRATO <- framecens$LABEL
    framecens$WEIGHTS <- 1
    # framecens$FPC <- 1
  }
  # Bethel multivariate allocation----------------------------------------
  strata <- strata[order(strata$STRATO),]
  framesamp <- framesamp[order(framesamp$STRATO),]
  newstratatot <- NULL
  for (j in (1:length(unique(strata$DOM1)))) {
    SOLUZ <- as.numeric(bethel(strata[strata$DOM1==j,],errors[j,],minnumstrat=minnumstrat))
    newstrata <- cbind(strata[strata$DOM1==j,],SOLUZ)
    newstratatot <- rbind(newstratatot,newstrata)
  }
  expected_CV(newstratatot)
  strata <- newstratatot
  
  # Sample selection -----------------------------------------------------
  
  if (sampling_method=="srs") samp <- selectSample(framesamp,strata) 
  if (sampling_method=="systematic") samp <- selectSampleSystematic(framesamp,strata) 
  if (sampling_method=="spatial") samp <- selectSampleSpatial(framesamp,strata,coord_names=c("LON","LAT")) 
  samp$FPC <- NULL
  
  # Output preparation -----------------------------------------------------
  
  samptot <- rbind(samp,framecens)
  outstrata <- rbind(strata,cens)
  framecens$WEIGHTS <- NULL
  framecens$FPC <- NULL
  # frametot <- rbind(framesamp,framecens)
  out <- list(sample=samptot,
              strata=strata,
              cens=cens,
              framesamp=framesamp,
              framecens=framecens)
  return(out)
}