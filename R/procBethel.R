procBethel <- function (framesamp, 
                        framecens, 
                        errors, 
                        sampling_method = c("srs","systematic", "spatial"), 
                        minnumstrat = 2) 
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
    cens <- buildStrataDF(framecens, progress = F)
    cens$CENS <- 1
    cens$SOLUZ <- cens$N
    colnames(framecens) <- toupper(colnames(framecens))
    cens$STRATUM <- as.numeric(as.factor(cens$STRATO))
    st1 <- "framecens$STRATUM <- paste(framecens$DOMAINVALUE,framecens$X"
    for (i in c(1:nvarX)) {
      if (i < nvarX) 
        st1 <- paste0(st1, i, ",framecens$X")
      if (i == nvarX) 
        st1 <- paste0(st1, i, ",sep='*')")
    }
    eval(parse(text = st1))
    framecens$LABEL <- as.numeric(as.factor(framecens$STRATUM))
    framecens$STRATO <- framecens$LABEL
    framecens$WEIGHTS <- 1
  }
  strata <- strata[order(strata$STRATO), ]
  framesamp <- framesamp[order(framesamp$STRATO), ]
  newstratatot <- NULL
  for (j in (1:length(unique(strata$DOM1)))) {
    SOLUZ <- as.numeric(bethel(strata[strata$DOM1 == j, 
    ], errors[j, ], minnumstrat = minnumstrat))
    newstrata <- cbind(strata[strata$DOM1 == j, ], SOLUZ)
    newstratatot <- rbind(newstratatot, newstrata)
  }
  expected_CV(newstratatot)
  strata <- newstratatot
  if (sampling_method == "srs") 
    samp <- selectSample(framesamp, strata)
  if (sampling_method == "systematic") 
    samp <- selectSampleSystematic(framesamp, strata)
  if (sampling_method == "spatial") 
    samp <- selectSampleSpatial(framesamp, strata, coord_names = c("LON","LAT"))
  samp$FPC <- NULL
  if (is.null(framecens)) {
    samptot <- samp
    outstrata <- strata
    out <- list(sample = samptot, strata = strata, framesamp = framesamp)
  }
  if (!is.null(framecens)) {
    samptot <- rbind(samp, framecens)
    outstrata <- rbind(strata, cens)
    framecens$WEIGHTS <- NULL
    framecens$FPC <- NULL
    out <- list(sample = samptot, strata = strata, cens = cens, 
                framesamp = framesamp, framecens = framecens)
  }
  return(out)
}
