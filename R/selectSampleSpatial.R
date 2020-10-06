#-----------------------------------------------------------
# Function to select a sample using the lpm2_kdtree function
#  from the SamplingBigData package (Till?-Grafstrom)
#  in order to spatially distribute selected points
#  for any stratum
#-----------------------------------------------------------

selectSampleSpatial <- function(frame, outstrata, coord_names)  {
  if (!(coord_names[1] %in% colnames(frame) )) stop("Coordinate names not in frame")
  if (!(coord_names[2] %in% colnames(frame) )) stop("Coordinate names not in frame")
  samp <- NULL
  for (i in unique(frame$DOMAINVALUE)) {
    # cat("\n Domain",i)
    strata <- outstrata[outstrata$DOM1 == i,]
    for (j in as.numeric(unique(strata$STRATO))) {
      strato <- frame[frame$DOMAINVALUE == i  & frame$LABEL == j,]
      st <- paste0("X <- cbind(strato$",coord_names[1],",strato$",coord_names[2],")")
      eval(parse(text=st))
      size <- round(strata$SOLUZ[strata$DOM1 == i & strata$STRATO == j])
      pop <- strata$N[strata$DOM1 == i & strata$STRATO == j]
      p <- rep( size / pop, nrow(X))
      s = lpm2_kdtree(p,X)
      s <- strato[s,]
      s$WEIGHTS <- nrow(strato) / nrow(s)
      samp <- rbind(samp,s)
    }
  }
  cat("\n*** Sample has been drawn successfully ***")
  cat("\n", nrow(samp), " units have been selected from ", 
      nrow(outstrata), " strata\n")
  return(samp)
}
