\name{selectSampleSpatial}
\Rdversion{1.5-2}
\alias{selectSampleSpatial}
\title{
Selection of geo-referenced points from the frame  
}
\description{
Once optimal stratification has been obtained, and a new frame has been 
built by assigning to the units of the old one the new strata labels, 
it is possible to select a stratified sample from the frame. 
If geographical coordinates are available in the frame, in order to obtain spatially distributed selected points
this function makes use of the 'lpm2_kdtree' function from the' SamplingBigData' package (Lisic-Grafstrom).

}
\usage{
selectSampleSpatial(frame, outstrata, coord_names)
}
\arguments{
  \item{frame}{
  This is the (mandatory) dataframe containing the sampling frame, as it has been modified 
  by the execution of the "updateFrame" function. 
}
  \item{outstrata}{
  This is the (mandatory) dataframe containing the information related to resulting stratification obtained by the execution of "optimStrata" function. It should coincide with 'solution$aggr_strata'.
}
  \item{coord_names}{
  Indicates with which names the coordinates are indicated in the frame.
}
}

\value{
A dataframe containing the selected sample
}

\author{
Giulio Barcaroli
}
 
\examples{
\dontrun{
#############################
# Example of "spatial" method
#############################
library(sp)
# locations (155 observed points)
data("meuse")
# grid of points (3103)
data("meuse.grid")
meuse.grid$id <- c(1:nrow(meuse.grid))
coordinates(meuse)<-c("x","y")
coordinates(meuse.grid)<-c("x","y")

## Kriging model
library(automap)
kriging_lead = autoKrige(log(lead) ~ dist, meuse, meuse.grid)
plot(kriging_lead,sp.layout = NULL, justPosition = TRUE)
kriging_zinc = autoKrige(log(zinc) ~ dist, meuse, meuse.grid)
plot(kriging_zinc, sp.layout = list(pts = list("sp.points", meuse)))
r2_lead <- 1 - kriging_lead$sserr/sum((meuse$lead-mean(meuse$lead))^2)
r2_lead
r2_zinc <- 1 - kriging_zinc$sserr/sum((meuse$zinc-mean(meuse$zinc))^2)
r2_zinc
df <- NULL
df$id <- meuse.grid$id
df$lead.pred <- kriging_lead$krige_output@data$var1.pred
df$lead.var <- kriging_lead$krige_output@data$var1.var
df$zinc.pred <- kriging_zinc$krige_output@data$var1.pred
df$zinc.var <- kriging_zinc$krige_output@data$var1.var
df$lon <- meuse.grid$x
df$lat <- meuse.grid$y
df$dom1 <- 1
df <- as.data.frame(df)
head(df)

## Optimization 
library(SamplingStrata)
frame <- buildFrameSpatial(df=df,
                           id="id",
                           X=c("lead.pred","zinc.pred"),
                           Y=c("lead.pred","zinc.pred"),
                           variance=c("lead.var","zinc.var"),
                           lon="lon",
                           lat="lat",
                           domainvalue = "dom1")
cv <- as.data.frame(list(DOM=rep("DOM1",1),
                         CV1=rep(0.01,1),
                         CV2=rep(0.01,1),
                         domainvalue=c(1:1) ))
set.seed(1234)
solution <- optimStrata (
  method = "spatial",
  errors=cv, 
  framesamp=frame,
  iter = 15,
  pops = 10,
  nStrata = 5,
  fitting = c(r2_lead,r2_zinc),
  range = c(kriging_lead$var_model$range[2],kriging_zinc$var_model$range[2]),
  kappa=1,
  writeFiles = FALSE,
  showPlot = TRUE,
  parallel = FALSE)
framenew <- solution$framenew
outstrata <- solution$aggr_strata

# Sample selection
samp <- selectSampleSpatial(framenew,outstrata,coord_names=c("LON","LAT"))
table(samp$STRATO)

# Plot
library(sf)
samp_sf <- st_as_sf(samp, coords = c("LON", "LAT"))
plot(samp_sf["STRATO"])

}
}
\keyword{ survey }

