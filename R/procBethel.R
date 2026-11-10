# ==============================================================================
# File:      procBethel.R
# Functions: procBethel() [exported]
#
# Purpose:
#   Complete procedure for given (not optimized) strata: the strata are the
#   cross-classification of domain and X variables of framesamp; it computes
#   the Bethel allocation, adds the take-all units in framecens, labels the
#   frames and selects the sample with the chosen method ("srs", "systematic"
#   or "spatial"), returning sample, strata, take-all strata and labelled
#   frames.
# ==============================================================================

#' @title Procedure to apply Bethel algorithm and select a sample from given strata
#'
#' @description
#' This function allows to execute a complete procedure from the Bethel optimal allocation to the selection of a sample, without having to optimize the strata, that are supposed to be given and fixed.
#'
#' @param framesamp Dataframe containing sampling frame units.
#' @param framecens Dataframe containing frame units that must be selected.
#' @param errors Dataframe containing the precision levels expressed in terms of maximum expected value of the Coefficients of Variation related to target variables of the survey.
#' @param sampling_method Parameter for choosing the selection method: "srs", "systematic" and "spatial".
#' @param minnumstrat Indicates the minimum number of units that must be allocated in each stratum. Default is 2.
#'
#' @return
#' List containing
#' (1) the selected sample,
#' (2) the strata with allocated sampling units,
#' (3) the take-all strata
#' (4) the sampling frame with the labels linking to (2)
#' (5) the frame with take-all units, with the labels linking to (3),
#' (6) the expected CVs of the target variables in each domain (element
#' \code{expected_cv}).
#' Elements (3) and (5) are present only when \code{framecens} is given.
#'
#' @author
#' Giulio Barcaroli
#'
#' @examples
#' \dontrun{
#' data("swissmunicipalities")
#' swissmun <- swissmunicipalities[swissmunicipalities$REG < 4,
#'                                 c("REG","COM","Nom","HApoly",
#'                                   "Surfacesbois","Surfacescult",
#'                                   "Airbat","POPTOT")]
#' ndom <- length(unique(swissmun$REG))
#' cv <- as.data.frame(list(DOM=rep("DOM1",ndom),
#'                          CV1=rep(0.10,ndom),
#'                          CV2=rep(0.10,ndom),
#'                          domainvalue=c(1:ndom) ))
#' cv
#' swissmun$HApoly.cat <- var.bin(swissmun$HApoly,15)
#' swissmun$POPTOT.cat <- var.bin(swissmun$POPTOT,15)
#' frame <- buildFrameDF(df = swissmun,
#'                          id = "COM",
#'                          X = c("POPTOT.cat","HApoly.cat"),
#'                          Y = c("Airbat","Surfacesbois"),
#'                          domainvalue = "REG")
#' summary(frame)
#' #----Selection of units to be censused from the frame
#' ind_framecens <- which(frame$X1 > 9)
#' framecens <- frame[ind_framecens,]
#' #----Selection of units to be sampled from the frame
#' # (complement to the previous)
#' framesamp <- frame[-ind_framecens,]
#'
#' a <- procBethel(framesamp,framecens,errors=cv,sampling_method="srs",minnumstrat=2)
#' head(a$sample)
#' expected_CV(a$strata)
#' }
#'
#' @keywords survey
#' @export
procBethel <- function (framesamp, 
                        framecens, 
                        errors, 
                        sampling_method = c("srs","systematic", "spatial"), 
                        minnumstrat = 2) 
{
  sampling_method <- match.arg(sampling_method)
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
    cens <- buildStrataDF(framecens, progress = F)
    colnames(cens) <- toupper(colnames(cens))
    cens$CENS <- 1
    cens$SOLUZ <- cens$N
    # take-all strata: same "domain*X1*...*Xk" key as the sampling strata,
    # numeric labels following those of the sampling strata (no overlap)
    xcols <- paste0("X", 1:nvarX)
    framecens$STRATUM <- do.call(paste, c(list(framecens$DOMAINVALUE), framecens[xcols], sep = "*"))
    cens$STRATUM <- do.call(paste, c(list(cens$DOM1), cens[xcols], sep = "*"))
    censkeys <- sort(unique(framecens$STRATUM))
    framecens$LABEL <- max(strata$STRATO) + match(framecens$STRATUM, censkeys)
    framecens$STRATO <- framecens$LABEL
    framecens$WEIGHTS <- 1
    cens$STRATO <- max(strata$STRATO) + match(cens$STRATUM, censkeys)
  }
  strata <- strata[order(strata$STRATO), ]
  framesamp <- framesamp[order(framesamp$STRATO), ]
  newstratatot <- NULL
  colnames(errors) <- toupper(colnames(errors))
  for (j in sort(unique(strata$DOM1))) {
    # precision constraints of domain j (not simply the j-th row of errors)
    SOLUZ <- as.numeric(bethel(strata[strata$DOM1 == j, 
    ], errors[errors$DOMAINVALUE == j, ], minnumstrat = minnumstrat))
    newstrata <- cbind(strata[strata$DOM1 == j, ], SOLUZ)
    newstratatot <- rbind(newstratatot, newstrata)
  }
  # expected CVs of the allocation (previously computed and discarded)
  expcv <- expected_CV(newstratatot)
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
    out <- list(sample = samptot, strata = strata, framesamp = framesamp,
                expected_cv = expcv)
  }
  if (!is.null(framecens)) {
    samptot <- rbind(samp, framecens)
    outstrata <- rbind(strata, cens)
    framecens$WEIGHTS <- NULL
    framecens$FPC <- NULL
    out <- list(sample = samptot, strata = strata, cens = cens, 
                framesamp = framesamp, framecens = framecens,
                expected_cv = expcv)
  }
  return(out)
}
