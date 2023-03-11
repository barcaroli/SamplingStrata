procBethel <- function (framesamp, framecens, errors, sampling_method = c("srs", 
                                                                          "systematic", "spatial"), minnumstrat = 2) 
{
  colnames(framesamp) <- toupper(colnames(framesamp))
  nvarX <- length(grep("X", colnames(framesamp)))
  nvarY <- length(grep("Y", colnames(framesamp)))
  st1 <- "framesamp$STRATUM <- paste(framesamp$DOMAINVALUE,framesamp$X"
  for (i in c(1:nvarX)) {
    if (i < nvarX) 
      st1 <- paste0(st1, i, ",framesamp$X")
    if (i == nvarX) 
      st1 <- paste0(st1, i, ",sep='*')")
  }
  eval(parse(text = st1))
  framesamp$STRATO <- as.numeric(as.factor(framesamp$STRATUM))
  framesamp$LABEL <- framesamp$STRATO
  strata <- buildStrataDF(framesamp, progress = F)
  colnames(strata) <- toupper(colnames(strata))
  st2 <- "strata$STRATUM <- paste(strata$DOM1,strata$X"
  for (i in c(1:nvarX)) {
    if (i < nvarX) 
      st2 <- paste0(st2, i, ",strata$X")
    if (i == nvarX) 
      st2 <- paste0(st2, i, ",sep='*')")
  }
  eval(parse(text = st2))
  strata$STRATO <- as.numeric(as.factor(strata$STRATUM))
  if (!is.null(framecens)) {
    colnames(framecens) <- toupper(colnames(framecens))
    nvarX <- length(grep("X", colnames(framecens)))
    nvarY <- length(grep("Y", colnames(framecens)))
    st1 <- "framecens$STRATUM <- paste(framecens$DOMAINVALUE,framecens$X"
    for (i in c(1:nvarX)) {
      if (i < nvarX) 
        st1 <- paste0(st1, i, ",framecens$X")
      if (i == nvarX) 
        st1 <- paste0(st1, i, ",sep='*')")
    }
    eval(parse(text = st1))
    framecens$STRATO <- as.numeric(as.factor(framecens$STRATUM))
    framecens$LABEL <- framecens$STRATO
    cens <- buildStrataDF(framecens, progress = F)
    colnames(cens) <- toupper(colnames(cens))
    st2 <- "cens$STRATUM <- paste(cens$DOM1,cens$X"
    for (i in c(1:nvarX)) {
      if (i < nvarX) 
        st2 <- paste0(st2, i, ",cens$X")
      if (i == nvarX) 
        st2 <- paste0(st2, i, ",sep='*')")
    }
    eval(parse(text = st2))
    cens$STRATO <- as.numeric(as.factor(cens$STRATUM))
    cens$CENS <- 1
    framecens$WEIGHTS <- 1
  }
  strata <- strata[order(strata$STRATO), ]
  framesamp <- framesamp[order(framesamp$STRATO), ]
  newstratatot <- NULL
  if (!is.null(framecens)) stratatot <- rbind(strata,cens)
  if (is.null(framecens)) stratatot <-strata
  for (j in (1:length(unique(stratatot$DOM1)))) {
    SOLUZ <- as.numeric(bethel(stratatot[stratatot$DOM1 == j, 
    ], errors[j, ], minnumstrat = minnumstrat))
    newstrata <- cbind(stratatot[stratatot$DOM1 == j, ], SOLUZ)
    newstratatot <- rbind(newstratatot, newstrata)
  }
  expected_CV(newstratatot)
  strata <- newstratatot[newstratatot$CENS == 0,]
  if (sampling_method == "srs") 
    samp <- selectSample(framesamp, strata)
  if (sampling_method == "systematic") 
    samp <- selectSampleSystematic(framesamp, strata)
  if (sampling_method == "spatial") 
    samp <- selectSampleSpatial(framesamp, strata, coord_names = c("LON","LAT"))
  samp$FPC <- NULL
  framecens$WEIGHTS <- 1
  samptot <- rbind(samp, framecens)
  cens$SOLUZ <- cens$N
  outstrata <- rbind(strata, cens)
  framecens$FPC <- NULL
  framecens$WEIGHTS <- NULL
  out <- list(sample = samptot, strata = outstrata, cens = cens, 
              framesamp = framesamp, framecens = framecens)
  return(out)
}
