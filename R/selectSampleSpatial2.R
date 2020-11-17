selectSampleSpatial2 <- function(frame, outstrata, coord_names)  {
  if (!(coord_names[1] %in% colnames(frame) )) stop("Coordinate names not in frame")
  if (!(coord_names[2] %in% colnames(frame) )) stop("Coordinate names not in frame")
  colnames(frame) <- toupper(colnames(frame))
  X <- cbind(frame$LON,frame$LAT)
  p <- rep( outstrata$SOLUZ / outstrata$N, outstrata$N)
  s = lpm2_kdtree(p,X)
  samp <- frame[s,]
  samp$WEIGHTS <- 1 / p[s]
  cat("\n*** Sample has been drawn successfully ***")
  cat("\n", nrow(samp), " units have been selected from ", 
      nrow(outstrata), " strata\n")
  return(samp)
}
