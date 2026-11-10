# ==============================================================================
# File:      selectSampleSpatial.R
# Functions: selectSampleSpatial() [exported]
#
# Purpose:
#   Selects a spatially balanced sample of geo-referenced units from the frame
#   labelled with the optimized strata, using the local pivotal method
#   (lpm2_kdtree() from SamplingBigData) within each stratum, with equal
#   inclusion probabilities SOLUZ / N; returns the sample with weights.
# ==============================================================================

#-----------------------------------------------------------
# Function to select a sample using the lpm2_kdtree function
#  from the SamplingBigData package (Tille-Grafstrom)
#  in order to spatially distribute selected points
#  for any stratum
#-----------------------------------------------------------

#' @title Selection of geo-referenced points from the frame
#'
#' @description
#' Once optimal stratification has been obtained, and a new frame has been
#' built by assigning to the units of the old one the new strata labels,
#' it is possible to select a stratified sample from the frame.
#' If geographical coordinates are available in the frame, in order to obtain spatially distributed selected points
#' this function makes use of the 'lpm2_kdtree' function from the' SamplingBigData' package (Lisic-Grafstrom).
#'
#' @param frame This is the (mandatory) dataframe containing the sampling frame, as it has been modified
#'   by the execution of the "updateFrame" function.
#' @param outstrata This is the (mandatory) dataframe containing the information related to resulting stratification obtained by the execution of "optimStrata" function. It should coincide with 'solution$aggr_strata'.
#' @param coord_names Indicates with which names the coordinates are indicated in the frame.
#'
#' @return
#' A dataframe containing the selected sample
#'
#' @author
#' Giulio Barcaroli
#'
#' @examples
#' \dontrun{
#' #############################
#' # Example of "spatial" method
#' #############################
#' library(sp)
#' # locations (155 observed points)
#' data("meuse")
#' # grid of points (3103)
#' data("meuse.grid")
#' meuse.grid$id <- c(1:nrow(meuse.grid))
#' coordinates(meuse)<-c("x","y")
#' coordinates(meuse.grid)<-c("x","y")
#'
#' ## Kriging model
#' library(automap)
#' kriging_lead = autoKrige(log(lead) ~ dist, meuse, meuse.grid)
#' plot(kriging_lead,sp.layout = NULL, justPosition = TRUE)
#' kriging_zinc = autoKrige(log(zinc) ~ dist, meuse, meuse.grid)
#' plot(kriging_zinc, sp.layout = list(pts = list("sp.points", meuse)))
#' r2_lead <- 1 - kriging_lead$sserr/sum((meuse$lead-mean(meuse$lead))^2)
#' r2_lead
#' r2_zinc <- 1 - kriging_zinc$sserr/sum((meuse$zinc-mean(meuse$zinc))^2)
#' r2_zinc
#' df <- NULL
#' df$id <- meuse.grid$id
#' df$lead.pred <- kriging_lead$krige_output@data$var1.pred
#' df$lead.var <- kriging_lead$krige_output@data$var1.var
#' df$zinc.pred <- kriging_zinc$krige_output@data$var1.pred
#' df$zinc.var <- kriging_zinc$krige_output@data$var1.var
#' df$lon <- meuse.grid$x
#' df$lat <- meuse.grid$y
#' df$dom1 <- 1
#' df <- as.data.frame(df)
#' head(df)
#'
#' ## Optimization
#' library(SamplingStrata)
#' frame <- buildFrameSpatial(df=df,
#'                            id="id",
#'                            X=c("lead.pred","zinc.pred"),
#'                            Y=c("lead.pred","zinc.pred"),
#'                            variance=c("lead.var","zinc.var"),
#'                            lon="lon",
#'                            lat="lat",
#'                            domainvalue = "dom1")
#' cv <- as.data.frame(list(DOM=rep("DOM1",1),
#'                          CV1=rep(0.01,1),
#'                          CV2=rep(0.01,1),
#'                          domainvalue=c(1:1) ))
#' set.seed(1234)
#' solution <- optimStrata (
#'   method = "spatial",
#'   errors=cv,
#'   framesamp=frame,
#'   iter = 15,
#'   pops = 10,
#'   nStrata = 5,
#'   fitting = c(r2_lead,r2_zinc),
#'   range = c(kriging_lead$var_model$range[2],kriging_zinc$var_model$range[2]),
#'   kappa=1,
#'   writeFiles = FALSE,
#'   showPlot = TRUE,
#'   parallel = FALSE)
#' framenew <- solution$framenew
#' outstrata <- solution$aggr_strata
#'
#' # Sample selection
#' samp <- selectSampleSpatial(framenew,outstrata,coord_names=c("LON","LAT"))
#' table(samp$STRATO)
#'
#' # Plot
#' library(sf)
#' samp_sf <- st_as_sf(samp, coords = c("LON", "LAT"))
#' plot(samp_sf["STRATO"])
#'
#' }
#'
#' @keywords survey
#' @export
selectSampleSpatial <- function(frame, outstrata, coord_names)  {
  if (!(coord_names[1] %in% colnames(frame) )) stop("Coordinate names not in frame")
  if (!(coord_names[2] %in% colnames(frame) )) stop("Coordinate names not in frame")
  X <- NULL
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
