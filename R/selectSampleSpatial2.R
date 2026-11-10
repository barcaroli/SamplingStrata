# ==============================================================================
# File:      selectSampleSpatial2.R
# Functions: selectSampleSpatial2() [internal, not exported]
#
# Purpose:
#   Earlier variant of selectSampleSpatial(): local pivotal method
#   (lpm2_kdtree()) on the whole frame, with coordinates taken from the columns
#   indicated in coord_names.
# ==============================================================================

selectSampleSpatial2 <- function(frame, outstrata, coord_names)  {
  if (!(coord_names[1] %in% colnames(frame) )) stop("Coordinate names not in frame")
  if (!(coord_names[2] %in% colnames(frame) )) stop("Coordinate names not in frame")
  # coordinates taken from the columns indicated in coord_names (they were
  # fixed to LON and LAT)
  X <- cbind(frame[[coord_names[1]]], frame[[coord_names[2]]])
  colnames(frame) <- toupper(colnames(frame))
  colnames(outstrata) <- toupper(colnames(outstrata))
  # inclusion probability of each unit = sampling rate of its stratum, matched
  # on (domain, stratum) instead of assuming frame and outstrata in the same order
  key_frame <- paste(frame$DOMAINVALUE, frame$LABEL)
  key_strata <- paste(outstrata$DOM1, outstrata$STRATO)
  p <- (outstrata$SOLUZ / outstrata$N)[match(key_frame, key_strata)]
  if (anyNA(p)) stop("Some units of frame belong to strata not in outstrata")
  s = lpm2_kdtree(p,X)
  samp <- frame[s,]
  samp$WEIGHTS <- 1 / p[s]
  cat("\n*** Sample has been drawn successfully ***")
  cat("\n", nrow(samp), " units have been selected from ", 
      nrow(outstrata), " strata\n")
  return(samp)
}
